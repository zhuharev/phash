package main

// Go provides a `flag` package supporting basic
// command-line flag parsing. We'll use this package to
// implement our example command-line program.
import "flag"
import "os"
import "log"
import "image"
import "image/jpeg"
import "image/png"
import "github.com/zhuharev/phash"

func init() {
	// damn important or else At(), Bounds() functions will
	// caused memory pointer error!!
	image.RegisterFormat("jpeg", "jpeg", jpeg.Decode, jpeg.DecodeConfig)
	image.RegisterFormat("png", "png", png.Decode, png.DecodeConfig)
}

func main() {
	image1Ptr := flag.String("i1", "", "input image 1")
	image2Ptr := flag.String("i2", "", "input image 2")
	image1OutPtr := flag.String("o1", "", "output image 1")
	image2OutPtr := flag.String("o2", "", "output image 2")

	flag.Parse()

	if *image1Ptr == "" {
		log.Fatal("Need to provide i1 flag for path to image 1")
	}

	if *image2Ptr == "" {
		log.Fatal("Need to provide i2 flag for path to image 2")
	}

	image1File, err := os.Open(*image1Ptr)
	if err != nil {
		log.Fatalf("Failed to open image %s, %v", *image1Ptr, err)
	}
	defer image1File.Close()

	image2File, err := os.Open(*image2Ptr)
	if err != nil {
		log.Fatalf("Failed to open image %s, %v", *image2Ptr, err)
	}
	defer image2File.Close()

	im1, _, err := image.Decode(image1File)
	if err != nil {
		log.Fatalf("Failed to decode image %v", err)
	}

	im2, _, err := image.Decode(image2File)
	if err != nil {
		log.Fatalf("Failed to decode image %v", err)
	}

	hash1 := phash.GreyscaleDctMatrix(im1)
	hash2 := phash.GreyscaleDctMatrix(im2)

	distance := phash.HammingDistance(hash1, hash2)

	log.Printf("Image 1 hash: %x\n", hash1)
	log.Printf("Image 2 hash: %x\n", hash2)
	log.Printf("Hamming Distance: %d", distance)

	if *image1OutPtr != "" {
		im1 := phash.HashToImage(hash1)
		out, _ := os.Create(*image1OutPtr)
		var opt jpeg.Options
		opt.Quality = 100
		_ = jpeg.Encode(out, im1, &opt) // put quality to 80%
	}

	if *image2OutPtr != "" {
		im2 := phash.HashToImage(hash2)
		out, _ := os.Create(*image2OutPtr)
		var opt jpeg.Options
		opt.Quality = 100
		_ = jpeg.Encode(out, im2, &opt) // put quality to 80%
	}
}
