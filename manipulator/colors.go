package manipulator

import (
	"fmt"
	"image/color"
)

//ColorToFloat64 converts a color.Color to float64.
func ColorToFloat64(c color.Color) float64 {
	r, g, b, _ := c.RGBA()
	return float64((r>>8)<<16) + float64((g>>8)<<8) + float64(b>>8)
}

func Float64ToColor(src float64) color.Color {
	c := int64(src)
	return color.RGBA{uint8(c >> 16), uint8(c >> 8), uint8(c), 0xff}
}

//ColorToHexString converts an RGB triple to an Hex string.
func ColorToHexString(c color.Color) string {
	r, g, b, _ := c.RGBA()
	return fmt.Sprintf("0x%02X%02X%02X", uint8(r>>8), uint8(g>>8), uint8(b>>8))
}
