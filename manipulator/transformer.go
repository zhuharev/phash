package manipulator

import (
	"github.com/gonum/matrix/mat64"
	"image"
	"image/color"
	"image/draw"
	"math"
)

func CropMatrix(src mat64.Matrix, x0, y0, x1, y1 int) (*mat64.Dense, error) {
	rows := y1 - y0 + 1
	cols := x1 - x0 + 1

	mtx := make([]float64, rows*cols)
	for x := x0; x <= x1; x++ {
		for y := y0; y <= y1; y++ {
			mtx[(y-y0)*cols+(x-x0)] = src.At(x, y)
		}
	}
	return mat64.NewDense(rows, cols, mtx), nil
}

func GrayImageToMatrix(src image.Gray) (*mat64.Dense, error) {
	bounds := src.Bounds()
	rows := bounds.Max.Y
	cols := bounds.Max.X

	mtx := make([]float64, rows*cols)
	for x := 0; x < cols; x++ {
		for y := 0; y < rows; y++ {
			_, _, b, _ := src.At(x, y).RGBA()
			mtx[y*cols+x] = float64(b)
		}
	}
	return mat64.NewDense(rows, cols, mtx), nil
}

func ImageToMatrix(src image.Image) (*mat64.Dense, error) {
	bounds := src.Bounds()
	rows := bounds.Max.Y
	cols := bounds.Max.X

	mtx := make([]float64, rows*cols)
	for x := 0; x < cols; x++ {
		for y := 0; y < rows; y++ {
			mtx[y*cols+x] = ColorToFloat64(src.At(x, y))
		}
	}
	return mat64.NewDense(rows, cols, mtx), nil
}

//ImageToGrayscale returns the greyscale of an image
func ImageToGrayscale(src image.Image) image.Gray {
	// Create a new grayscale image
	bounds := src.Bounds()
	gray := image.NewGray(bounds)
	for x := 0; x < bounds.Max.X; x++ {
		for y := 0; y < bounds.Max.Y; y++ {
			gray.Set(x, y, color.GrayModel.Convert(src.At(x, y)))
		}
	}
	return *gray
}

func MatrixToImage(src mat64.Matrix) image.Image {
	width, height := src.Dims()
	img := image.NewRGBA(image.Rect(0, 0, width, height))
	for x := 0; x < width; x++ {
		for y := 0; y < height; y++ {
			img.Set(x, y, Float64ToColor(src.At(x, y)))
		}
	}
	return img
}

//CopyImage copies images into a draw.Image
func CopyImage(src image.Image) draw.Image {
	b := src.Bounds()
	copy := image.NewRGBA(b)

	draw.Draw(copy, copy.Bounds(), src, b.Min, draw.Src)

	return copy
}

const x = math.Pi / 180

func Rad(d float64) float64 {
	return d * x
}

func Deg(r float64) float64 {
	return r / x
}
