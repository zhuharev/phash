package phash

import (
	"fmt"
	"image"
	"image/color"
	// "github.com/hawx/img/greyscale"
	"github.com/disintegration/gift"
	"github.com/gonum/matrix/mat64"

	// Package image/[jpeg|fig|png] is not used explicitly in the code below,
	// but is imported for its initialization side-effect, which allows
	// image.Decode to understand [jpeg|gif|png] formatted images.
	_ "github.com/BurntSushi/graphics-go/graphics"
	"github.com/zhuharev/phash/manipulator"
	"github.com/zhuharev/phash/radon"
	_ "golang.org/x/image/bmp"
	_ "golang.org/x/image/tiff"
	// _ "github.com/kavu/go-phash"
	_ "image/gif"

	"github.com/nfnt/resize"
	_ "github.com/smartystreets/goconvey/convey"
	//"image/jpeg"
	_ "image/png"
	"math"
	//"os"
	"sort"
)

type ImageDigest struct {
	Radon       radon.ImageDigest
	Phash       uint64
	PhashMatrix uint64
}

func median(a []float64) float64 {
	b := make([]float64, len(a))
	copy(b, a)
	sort.Float64s(b)
	return b[int(len(b)/2)]
}

func coefficient(n int) float64 {
	if n == 0 {
		return 1.0 / math.Sqrt(2)
	}
	return 1.0
}

func dctPoint(img image.Gray, u, v, N, M int) float64 {
	sum := 0.0
	for i := 0; i < N; i++ {
		for j := 0; j < M; j++ {
			_, _, b, _ := img.At(i, j).RGBA()
			// sum += math.Cos( ( float64(2*i+1)/float64(2*N) ) * float64( u ) * math.Pi ) *
			//        math.Cos( ( float64(2*j+1)/float64(2*M) ) * float64( v ) * math.Pi ) *
			//        float64(b)

			sum += math.Cos(math.Pi/(float64(N))*(float64(i)+1/2)*float64(u)) *
				math.Cos(math.Pi/(float64(M))*(float64(j)+1/2)*float64(v)) *
				float64(b)
		}
	}
	return sum * ((coefficient(u) * coefficient(v)) / 4.0)
}

// GreyscaleDct Computes the Dct of a greyscale image
func GreyscaleDct(img image.Gray) uint64 {
	// func DctImageHashOne(img image.Image) ([][]float64) {
	R := img.Bounds()
	N := R.Dx() // width
	M := R.Dy() // height
	DCTMatrix := make([][]float64, N)
	for u := 0; u < N; u++ {
		DCTMatrix[u] = make([]float64, M)
		for v := 0; v < M; v++ {
			DCTMatrix[u][v] = dctPoint(img, u, v, N, M)
			// fmt.Println( "DCTMatrix[", u, "][", v, "] is ", DCTMatrix[u][v])
		}
	}

	total := 0.0
	for u := 0; u < N/2; u++ {
		for v := 0; v < M/2; v++ {
			total += DCTMatrix[u][v]
		}
	}
	total -= DCTMatrix[0][0]
	avg := total / float64(((N/2)*(M/2))-1)
	fmt.Println("got average ", avg)
	var hash uint64
	for u := 0; u < N/2; u++ {
		for v := 0; v < M/2; v++ {
			hash = hash * 2
			if DCTMatrix[u][v] > avg {
				hash++
			}
		}
	}

	return hash
}

func createDctMatrix(N, M int) (*mat64.Dense, error) {
	rows := N
	cols := M
	mtx := make([]float64, rows*cols)
	c1 := math.Sqrt(2.0 / float64(N))
	c0 := 1 / math.Sqrt(float64(M))
	for x := 0; x < cols; x++ {
		//mtx[x] = dctMatrixRow(N, M, x, c0, c1)

		//row := make([]float64, M)
		//row[0] = c0
		mtx[x] = c0
		for y := 1; y < rows; y++ {
			mtx[y*cols+x] = c1 * math.Cos((math.Pi/2.0/float64(N))*float64(y)*(2.0*float64(x)+1.0))
		}
	}

	return mat64.NewDense(rows, cols, mtx), nil
}

type FloatMatrix [][]float64

func NewFloatMatrix(rows int, columns int) FloatMatrix {
	r := make([][]float64, rows)
	for ii := 0; ii < rows; ii++ {
		r[ii] = make([]float64, columns)
	}
	return FloatMatrix(r)
}

func (m FloatMatrix) Rows() int {
	return len(m)
}

func (m FloatMatrix) Columns() int {
	if len(m) > 0 {
		return len(m[0])
	}
	return 0
}

func (m FloatMatrix) Transposed() FloatMatrix {
	t := NewFloatMatrix(m.Columns(), m.Rows())
	for ii := 0; ii < m.Rows(); ii++ {
		for jj := 0; jj < m.Columns(); jj++ {
			t[ii][jj] = m[jj][ii]
		}
	}
	return t
}

func (m FloatMatrix) SubMatrix(x int, y int, rows int, columns int) (FloatMatrix, error) {
	if x+rows > m.Rows() {
		return nil, fmt.Errorf("can't extract %d rows starting at %d, matrix has %d rows", rows, x, m.Rows())
	}
	if y+columns > m.Columns() {
		return nil, fmt.Errorf("can't extract %d columns starting at %d, matrix has %d columns", columns, y, m.Columns())
	}
	s := NewFloatMatrix(rows, columns)
	for ii := 0; ii < rows; ii++ {
		for jj := 0; jj < columns; jj++ {
			s[ii][jj] = m[ii+x][jj+y]
		}
	}
	return s, nil
}

func (m FloatMatrix) Multiply(n FloatMatrix) (FloatMatrix, error) {
	if m.Columns() != n.Rows() {
		return nil, fmt.Errorf("can't multiply matrix with %d columns by matrix with %d rows", m.Columns(), n.Rows())
	}
	p := NewFloatMatrix(m.Rows(), n.Columns())
	end := m.Columns()
	for ii := 0; ii < p.Rows(); ii++ {
		for jj := 0; jj < p.Columns(); jj++ {
			for kk := 0; kk < end; kk++ {
				p[ii][jj] += m[ii][kk] * n[kk][jj]
			}
		}
	}
	return p, nil
}

// UnrollX returns all values in a slice which contains
// all rows in increasing order.
func (m FloatMatrix) UnrollX() []float64 {
	values := make([]float64, 0, m.Rows()*m.Columns())
	for ii := 0; ii < m.Columns(); ii++ {
		for jj := 0; jj < m.Rows(); jj++ {
			values = append(values, m[jj][ii])
		}
	}
	return values
}

// UnrollX returns all values in a slice which contains
// all columns in increasing order.
func (m FloatMatrix) UnrollY() []float64 {
	values := make([]float64, 0, m.Rows()*m.Columns())
	for _, v := range m {
		values = append(values, v...)
	}
	return values
}

func NewDCTMatrix(order int) FloatMatrix {
	m := NewFloatMatrix(order, order)
	c0 := 1 / math.Sqrt(float64(order))
	c1 := math.Sqrt(2 / float64(order))
	for ii := 0; ii < order; ii++ {
		m[ii][0] = c0
		for jj := 1; jj < order; jj++ {
			m[ii][jj] = c1 * math.Cos((math.Pi/2/float64(order))*float64(jj)*(2*float64(ii)+1))
		}
	}
	return m
}

// GreyscaleDctMatrix Computes the Dct of a greyscale image using matrixes
func GreyscaleDctMatrix(im image.Image) uint64 {
	greyFilter := gift.Grayscale()

	convFilter := gift.Convolution([]float32{
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1,
	}, true, false, false, 0)
	resizeFilter := gift.Resize(32, 32, gift.LinearResampling)
	g := gift.New(greyFilter, convFilter, resizeFilter)
	dst := image.NewRGBA(g.Bounds(im.Bounds()))
	g.Draw(dst, im)
	/*
		out, _ := os.Create("/Users/danielriley/desktop/output.jpg")
		var opt jpeg.Options
		opt.Quality = 80
		_ = jpeg.Encode(out, dst, &opt) // put quality to 80%

		out1, _ := os.Create("/Users/danielriley/desktop/output1.jpg")
		opt.Quality = 80
		_ = jpeg.Encode(out1, im, &opt) // put quality to 80%
		return 0
	*/
	//width := im.Bounds().Max.X
	//height := im.Bounds().Max.Y
	m := make([][]float64, 32)
	for i := 0; i < 32; i++ {
		m[i] = make([]float64, 32)
		for j := 0; j < 32; j++ {
			_, _, b, _ := dst.At(i, j).RGBA()
			m[i][j] = float64(b)
		}
	}
	/*
		out, _ := os.Create("/Users/danielriley/desktop/output.jpg")
		var opt jpeg.Options
		opt.Quality = 80
		_ = jpeg.Encode(out, dst, &opt) // put quality to 80%

		out1, _ := os.Create("/Users/danielriley/desktop/output1.jpg")
		opt.Quality = 80
		_ = jpeg.Encode(out1, im, &opt) // put quality to 80%
	*/
	imMatrix := FloatMatrix(m)
	dctMatrix := NewDCTMatrix(32)

	// We can safely ignore errors here, since the sizes
	// always match.
	dct, _ := dctMatrix.Multiply(imMatrix)

	dct, _ = dct.Multiply(dctMatrix.Transposed())

	sub, _ := dct.SubMatrix(1, 1, 8, 8)

	values := sub.UnrollX()

	// We need the original values, so we must sort a copy
	cpy := make([]float64, len(values))
	copy(cpy, values)
	sort.Float64s(cpy)
	median := (cpy[64/2-1] + cpy[64/2]) / 2
	fmt.Println("got median ", median)
	bit := uint64(1)
	hash := uint64(0)
	for _, v := range values {
		if v > median {
			hash |= bit
		}
		bit <<= 1
	}

	fmt.Println("calculated hash ", hash)
	return hash

	/*
		imgMtx, err := manipulator.GrayImageToMatrix(img)
		if err != nil {
			panic(err)
		}
		dctMtx, err := createDctMatrix(img.Bounds().Max.X, img.Bounds().Max.Y)
		if err != nil {
			panic(err)
		}

		dctMtxTransp := dctMtx.T() // Transpose

		dctMtx.Mul(dctMtx, imgMtx)
		dctMtx.Mul(dctMtx, dctMtxTransp)

		dctImage, err := manipulator.CropMatrix(dctMtx, 0, 0, 7, 7)
		if err != nil {
			panic(err)
		}
		subsec := dctImage.RawMatrix().Data
		median := median(subsec)
		var one, hash uint64 = 1, 0
		for i := 0; i < len(subsec); i++ {
			current := subsec[i]
			if current > median {
				hash |= one
			}
			one = one << 1
		}
		return hash
	*/
}

func HashToImage(hash uint64) image.Image {
	im := image.NewRGBA(image.Rect(0, 0, 8, 8))

	bit := uint64(1)

	for i := 0; i < 64; i++ {
		if hash&bit != 0 {
			im.Set(i%8, i/8, color.RGBA{0, 0, 0, 255})
		} else {
			im.Set(i%8, i/8, color.RGBA{255, 255, 255, 255})
		}

		bit <<= 1
	}

	return im
}

//ComputeGreyscaleDct puts the result of GreyscaleDct in a digest
func (d *ImageDigest) ComputeGreyscaleDct() error {
	stamp := resize.Resize(32, 32, d.Radon.Image, resize.Bilinear)
	greyscaleStamp := manipulator.ImageToGrayscale(stamp)

	// greyscaleStamp := greyscale.Greyscale(stamp)
	d.Phash = GreyscaleDct(greyscaleStamp)

	return nil
}

//ComputeGreyscaleDctMatrix puts the result of GreyscaleDctMatrix in a digest
func (d *ImageDigest) ComputeGreyscaleDctMatrix() error {
	//stamp := resize.Resize(32, 32, d.Radon.Image, resize.Bilinear)
	//greyscaleStamp := manipulator.ImageToGrayscale(stamp)

	// greyscaleStamp := greyscale.Greyscale(stamp)
	d.PhashMatrix = GreyscaleDctMatrix(d.Radon.Image)

	return nil
}
