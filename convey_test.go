package phash_test

import (
	"fmt"
	// "github.com/azer-/phash"
	"code.google.com/p/graphics-go/graphics"
	. "github.com/smartystreets/goconvey/convey"
	"github.com/verisart/phash"
	"github.com/verisart/phash/manipulator"
	"github.com/verisart/phash/radon"
	"testing"
)

var testImages = map[string]*ImageBag{}

var Treshold = uint64(17)

func hammingDist(a, b *ImageBag) uint64 {
	dist := phash.HammingDistance(a.PhashMatrix, b.PhashMatrix)
	fmt.Printf("dist between %s and %s is %d", a.Filename, b.Filename, dist)
	return dist
}

func getImgBag(dir, filename string, angle Angle) *ImageBag {
	img, found := testImages[dir+filename]
	if !found {
		testImages[dir+filename] = new(ImageBag)
		testImages[dir+filename].Dir = dir
		testImages[dir+filename].Filename = filename
		testImages[dir+filename].Rotations = map[Angle]*ImageBag{}
		testImages[dir+filename].InitialiseFromFileInfo()
		img = testImages[dir+filename]
	}

	if angle != 0 {
		rotatedImage, found := img.Rotations[angle]
		if !found {
			draw := manipulator.CopyImage(img.Radon.Image)
			err := graphics.Rotate(draw, img.Radon.Image, &graphics.RotateOptions{Angle: float64(angle)})
			if err != nil {
				panic(err)
			}
			rImg := ImageBag{
				Dir:      dir,
				Filename: filename,
				Angle:    angle,
				ImageDigest: phash.ImageDigest{
					Radon: radon.ImageDigest{
						Image:  draw,
						Format: img.Radon.Format}}}
			img.Rotations[angle] = &rImg

			return img.Rotations[angle]
		}
		return rotatedImage
	}

	return img
}

func PairExecuteFor2ImagesList(imagesA []*ImageBag, imagesB []*ImageBag, f func(*ImageBag, *ImageBag)) {
	for _, imageA := range imagesA {
		for _, imageB := range imagesB {
			fmt.Println("checking ", imageA.Filename, imageB.Filename)
			f(imageA, imageB)
		}
	}
}

func TestAffineTransformedImagesMatch(t *testing.T) {

	Convey("Given some cats are loaded", t, func() {

		catBig := getImgBag(catsDir, "cat_big.jpg", 0)
		catMedium := getImgBag(catsDir, "cat_medium.jpg", 0)
		catSmall := getImgBag(catsDir, "cat_small.jpg", 0)
		/*
			Convey("When DCT Hashes are computed", func() {
				catBig.ImageDigest.ComputeGreyscaleDct()
				catMedium.ImageDigest.ComputeGreyscaleDct()
				catSmall.ImageDigest.ComputeGreyscaleDct()

				Convey("Then the hashes should not be zero", func() {
					So(catBig.Phash, ShouldNotEqual, 0)
					So(catMedium.Phash, ShouldNotEqual, 0)
					So(catSmall.Phash, ShouldNotEqual, 0)
				})

				Convey("And Hamming Distance should be under threshold ", func() {
					distance := phash.HammingDistance(catBig.Phash, catMedium.Phash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
					distance = phash.HammingDistance(catBig.Phash, catSmall.Phash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
					distance = phash.HammingDistance(catMedium.Phash, catSmall.Phash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
				})

			})
		*/
		Convey("When DCTMatrix Hashes are computed", func() {
			catBig.ImageDigest.ComputeGreyscaleDctMatrix()
			catMedium.ImageDigest.ComputeGreyscaleDctMatrix()
			catSmall.ImageDigest.ComputeGreyscaleDctMatrix()

			Convey("Then the hashes should not be zero", func() {
				So(catBig.PhashMatrix, ShouldNotEqual, 0)
				So(catMedium.PhashMatrix, ShouldNotEqual, 0)
				So(catSmall.PhashMatrix, ShouldNotEqual, 0)
			})

			Convey("Then Hamming Distance should be under threshold ", func() {
				distance := hammingDist(catBig, catMedium)
				So(distance, ShouldBeLessThanOrEqualTo, Treshold)
				distance = hammingDist(catBig, catSmall)
				So(distance, ShouldBeLessThanOrEqualTo, Treshold)
				distance = hammingDist(catMedium, catSmall)
				So(distance, ShouldBeLessThanOrEqualTo, Treshold)
			})

		})
		/*
			Convey("When CDCT Hashes are computed", func() {
				catBig.ComputeImageHashPhash(false)
				catMedium.ComputeImageHashPhash(false)
				catSmall.ComputeImageHashPhash(false)

				Convey("Then the hashes should not be zero", func() {
					So(catBig.CPhash, ShouldNotEqual, 0)
					So(catMedium.CPhash, ShouldNotEqual, 0)
					So(catSmall.CPhash, ShouldNotEqual, 0)
				})

				Convey("Then Hamming Distance should be under threshold ", func() {
					distance := phash.HammingDistance(catBig.CPhash, catMedium.CPhash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
					distance = phash.HammingDistance(catBig.CPhash, catSmall.CPhash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
					distance = phash.HammingDistance(catMedium.CPhash, catSmall.CPhash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
				})

			})
		*/
		SkipConvey("When Radon Hashes are computed", func() {
			catBig.ComputeImageHashRadon(false)
			catMedium.ComputeImageHashRadon(false)

			Convey("Then the hashes should not be zero", func() {
				So(len(catBig.ImageDigest.Radon.Digest.Coeffs), ShouldNotEqual, 0)
				So(catBig.ImageDigest.Radon.Digest.Coeffs, ShouldNotBeNil)
				So(len(catMedium.ImageDigest.Radon.Digest.Coeffs), ShouldNotEqual, 0)
				So(catMedium.ImageDigest.Radon.Digest.Coeffs, ShouldNotBeNil)
			})

			Convey("Then the Cross Correlation for threshold should be true", func() {
				crossCoorrelation, _ := radon.DiffByCrossCorr(catBig.ImageDigest.Radon.Digest.Coeffs, catMedium.ImageDigest.Radon.Digest.Coeffs, 0.01)
				So(crossCoorrelation, ShouldBeTrue)
			})

		})

	})

	Convey("Given some suns are loaded", t, func() {
		sunBig := getImgBag(catsDir, "sun_big.jpg", 0)
		sunSmall := getImgBag(catsDir, "sun_small.jpg", 0)

		Convey("When DCT Hashes are computed", func() {
			sunBig.ImageDigest.ComputeGreyscaleDct()
			sunSmall.ImageDigest.ComputeGreyscaleDct()

			Convey("Then the hashes should not be zero", func() {
				So(sunBig.Phash, ShouldNotEqual, 0)
				So(sunSmall.Phash, ShouldNotEqual, 0)
			})

			Convey("And Hamming Distance should be under threshold ", func() {
				distance := phash.HammingDistance(sunBig.Phash, sunSmall.Phash)
				So(distance, ShouldBeLessThanOrEqualTo, Treshold)
			})

		})

		Convey("When DCTMatrix Hashes are computed", func() {
			fmt.Println("checking sunBig")
			sunBig.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("checking sunSmall")
			sunSmall.ImageDigest.ComputeGreyscaleDctMatrix()

			Convey("Then the hashes should not be zero", func() {
				So(sunBig.PhashMatrix, ShouldNotEqual, 0)
				So(sunSmall.PhashMatrix, ShouldNotEqual, 0)
			})

			Convey("Then Hamming Distance should be under threshold ", func() {
				distance := hammingDist(sunBig, sunSmall)
				So(distance, ShouldBeLessThanOrEqualTo, Treshold)
				distance = hammingDist(sunBig, sunSmall)
				So(distance, ShouldBeLessThanOrEqualTo, Treshold)
			})

		})
		/*
			Convey("When CDCT Hashes are computed", func() {
				sunBig.ComputeImageHashPhash(false)
				sunSmall.ComputeImageHashPhash(false)

				Convey("Then the hashes should not be zero", func() {
					So(sunBig.CPhash, ShouldNotEqual, 0)
					So(sunSmall.CPhash, ShouldNotEqual, 0)
				})

				Convey("Then Hamming Distance should be under threshold ", func() {
					distance := phash.HammingDistance(sunBig.CPhash, sunSmall.CPhash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
					distance = phash.HammingDistance(sunBig.CPhash, sunSmall.CPhash)
					So(distance, ShouldBeLessThanOrEqualTo, Treshold)
				})

			})*/

		Convey("When Radon Hashes are computed", func() {
			sunBig.ComputeImageHashRadon(false)
			sunSmall.ComputeImageHashRadon(false)

			Convey("Then the hashes should not be zero", func() {
				So(len(sunBig.ImageDigest.Radon.Digest.Coeffs), ShouldNotEqual, 0)
				So(sunBig.ImageDigest.Radon.Digest.Coeffs, ShouldNotBeNil)
				So(len(sunSmall.ImageDigest.Radon.Digest.Coeffs), ShouldNotEqual, 0)
				So(sunSmall.ImageDigest.Radon.Digest.Coeffs, ShouldNotBeNil)
			})

			SkipConvey("Then the Cross Correlation for threshold should be true", func() {

				crossCoorrelation, _ := radon.DiffByCrossCorr(sunBig.ImageDigest.Radon.Digest.Coeffs, sunSmall.ImageDigest.Radon.Digest.Coeffs, -0.1)
				So(crossCoorrelation, ShouldBeTrue)
				crossCoorrelation, _ = radon.DiffByCrossCorr(sunBig.ImageDigest.Radon.Digest.Coeffs, sunSmall.ImageDigest.Radon.Digest.Coeffs, 0.0)
				So(crossCoorrelation, ShouldBeTrue)
			})

		})

	})

}

func TestDifferentImagesDoNotMatch(t *testing.T) {

	Convey("Given some cats, suns and girls are loaded", t, func() {

		catBig := getImgBag(catsDir, "cat_big.jpg", 0)
		catMedium := getImgBag(catsDir, "cat_medium.jpg", 0)
		catSmall := getImgBag(catsDir, "cat_small.jpg", 0)
		catSky := getImgBag(catsDir, "cat_sky.jpg", 0)
		catSmile := getImgBag(catsDir, "smiling.jpg", 0)
		cats := []*ImageBag{catBig, catMedium, catSmall, catSky, catSmile}

		lenaToshop := getImgBag(lenaDir, "lena_std.jpg", 0)
		lenaNkd := getImgBag(lenaDir, "l_hires.jpg", 0)
		lenaVintage := getImgBag(lenaDir, "lena.bmp", 0)
		girls := []*ImageBag{lenaToshop, lenaNkd, lenaVintage}

		sunBig := getImgBag(catsDir, "sun_big.jpg", 0)
		sunSmall := getImgBag(catsDir, "sun_small.jpg", 0)
		suns := []*ImageBag{sunBig, sunSmall}

		Convey("When DCT Hashes are computed", func() {

			catBig.ImageDigest.ComputeGreyscaleDct()
			catMedium.ImageDigest.ComputeGreyscaleDct()
			catSmall.ImageDigest.ComputeGreyscaleDct()
			catSky.ImageDigest.ComputeGreyscaleDct()
			catSmile.ImageDigest.ComputeGreyscaleDct()

			lenaToshop.ImageDigest.ComputeGreyscaleDct()
			lenaNkd.ImageDigest.ComputeGreyscaleDct()
			lenaVintage.ImageDigest.ComputeGreyscaleDct()

			sunBig.ImageDigest.ComputeGreyscaleDct()
			sunSmall.ImageDigest.ComputeGreyscaleDct()

			Convey("Then cats do not look like girls", func() {
				PairExecuteFor2ImagesList(cats, girls, func(cat, girl *ImageBag) {
					distance := phash.HammingDistance(cat.Phash, girl.Phash)
					So(distance, ShouldBeGreaterThan, Treshold)
				})
			})

			Convey("And suns do not look like cats", func() {
				PairExecuteFor2ImagesList(suns, cats, func(sun, cat *ImageBag) {
					distance := phash.HammingDistance(sun.Phash, cat.Phash)
					So(distance, ShouldBeGreaterThan, Treshold)
				})
			})

			Convey("And girls do not look like suns", func() {
				PairExecuteFor2ImagesList(girls, suns, func(girl, sun *ImageBag) {
					distance := phash.HammingDistance(sun.Phash, girl.Phash)
					So(distance, ShouldBeGreaterThan, Treshold)
				})
			})

		})

		Convey("When DCTMatrix Hashes are computed", func() {

			fmt.Println("catBig")
			catBig.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("catMedium")
			catMedium.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("catSmall")
			catSmall.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("catSky")
			catSky.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("catSmile")
			catSmile.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("lenaToshop")
			lenaToshop.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("lenankd")
			lenaNkd.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("lenavintage")
			lenaVintage.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("sunbig")
			sunBig.ImageDigest.ComputeGreyscaleDctMatrix()
			fmt.Println("sunsmall")
			sunSmall.ImageDigest.ComputeGreyscaleDctMatrix()

			Convey("Then cats do not look like girls", func() {
				PairExecuteFor2ImagesList(cats, girls, func(cat, girl *ImageBag) {
					distance := hammingDist(cat, girl)
					So(distance, ShouldBeGreaterThan, Treshold)
				})
			})

			Convey("And suns do not look like cats", func() {
				PairExecuteFor2ImagesList(suns, cats, func(sun, cat *ImageBag) {
					distance := hammingDist(sun, cat)
					So(distance, ShouldBeGreaterThan, Treshold)
				})
			})

			Convey("And girls do not look like suns", func() {
				PairExecuteFor2ImagesList(girls, suns, func(girl, sun *ImageBag) {
					distance := hammingDist(sun, girl)
					So(distance, ShouldBeGreaterThan, Treshold)
				})
			})

		})
		/*
			Convey("When CDCT Hashes are computed", func() {

				catBig.ComputeImageHashPhash(false)
				catMedium.ComputeImageHashPhash(false)
				catSmall.ComputeImageHashPhash(false)
				catSky.ComputeImageHashPhash(false)
				catSmile.ComputeImageHashPhash(false)

				// lenaToshop.ComputeImageHashPhash(false)
				lenaVintage.ComputeImageHashPhash(false)
				lenaNkd.ComputeImageHashPhash(false)
				girls := []*ImageBag{lenaVintage, lenaNkd}

				sunBig.ComputeImageHashPhash(false)
				sunSmall.ComputeImageHashPhash(false)

				Convey("Then cats do not look like girls", func() {
					PairExecuteFor2ImagesList(cats, girls, func(cat, girl *ImageBag) {
						distance := phash.HammingDistance(cat.CPhash, girl.CPhash)
						So(distance, ShouldBeGreaterThan, Treshold)
					})
				})

				Convey("And suns do not look like cats", func() {
					PairExecuteFor2ImagesList(suns, cats, func(sun, cat *ImageBag) {
						distance := phash.HammingDistance(sun.CPhash, cat.CPhash)
						So(distance, ShouldBeGreaterThan, Treshold)
					})
				})

				Convey("And girls do not look like suns", func() {
					PairExecuteFor2ImagesList(girls, suns, func(girl, sun *ImageBag) {
						distance := phash.HammingDistance(sun.CPhash, girl.CPhash)
						So(distance, ShouldBeGreaterThan, Treshold)
					})
				})

			})*/

		SkipConvey("When Radon Hashes are computed", func() {

			catBig.ComputeImageHashRadon(false)
			catMedium.ComputeImageHashRadon(false)
			catSmall.ComputeImageHashRadon(false)
			catSky.ComputeImageHashRadon(false)
			catSmile.ComputeImageHashRadon(false)

			lenaToshop.ComputeImageHashRadon(false)
			lenaNkd.ComputeImageHashRadon(false)
			lenaVintage.ComputeImageHashRadon(false)

			sunBig.ComputeImageHashRadon(false)
			sunSmall.ComputeImageHashRadon(false)

			Convey("Then cats do not look like girls", func() {
				PairExecuteFor2ImagesList(cats, girls, func(cat, girl *ImageBag) {
					crossCoorrelation, _ := radon.DiffByCrossCorr(cat.ImageDigest.Radon.Digest.Coeffs, girl.ImageDigest.Radon.Digest.Coeffs, 0.8)
					So(crossCoorrelation, ShouldBeFalse)
				})
			})

			Convey("And suns do not look like cats", func() {
				PairExecuteFor2ImagesList(suns, cats, func(sun, cat *ImageBag) {
					crossCoorrelation, _ := radon.DiffByCrossCorr(sun.ImageDigest.Radon.Digest.Coeffs, cat.ImageDigest.Radon.Digest.Coeffs, 0.8)
					So(crossCoorrelation, ShouldBeFalse)
				})
			})

			Convey("And girls do not look like suns", func() {
				PairExecuteFor2ImagesList(girls, suns, func(girl, sun *ImageBag) {
					crossCoorrelation, _ := radon.DiffByCrossCorr(girl.ImageDigest.Radon.Digest.Coeffs, sun.ImageDigest.Radon.Digest.Coeffs, 0.8)
					So(crossCoorrelation, ShouldBeFalse)
				})
			})

		})

	})

}
