// Copyright 2019 Hong-Ping Lo. All rights reserved.
// Use of this source code is governed by a BDS-style
// license that can be found in the LICENSE file.
package main

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"unsafe"

	"github.com/jblindsay/lidario"
	"gitlab.com/hongping1224/pointdensity/counter"
)

func TestFloat32Change(t *testing.T) {
	a := float32(100.0)
	fmt.Println(a)
	b := (*[4]byte)(unsafe.Pointer(&a))[:]
	fmt.Println(b)
	bits := binary.LittleEndian.Uint32(b)
	u := math.Float32frombits(bits)
	fmt.Println(u)
}

func TestSaveTiff(t *testing.T) {
	SaveTiff("./test", []int{70000}, 1, 1, 1, counter.Point{X: 0, Y: 0})
}

func TestPrintMap(t *testing.T) {
	numOFCPU = 4
	Run("testdata/sample.las", 0.1, "testdata/Output")
}

func TestOutputPath(t *testing.T) {
	inpath := "a\\b\\c\\123.las"

	empty := ""
	onlydir, _ := os.Getwd()
	onlydir += "\\"
	onlyfilename := "test.tiff"
	fullpath := onlydir + onlyfilename

	emptyFileAns := "a\\b\\c\\PointDensity_123"
	onlydirFileAns := filepath.Join(onlydir, "PointDensity_123")
	onlyfilenameFileAns := "a\\b\\c\\test"
	fullpathFileAns := strings.TrimSuffix(fullpath, ".tiff")
	emptyDirAns := "a\\b\\c\\PointDensity"
	onlydirDirAns := filepath.Join(onlydir, "PointDensity")
	onlyfilenameDirAns := "a\\b\\c\\PointDensity"
	fullpathDirAns := filepath.Join(onlydir, "PointDensity")
	//check filemode empty
	out := checkOutPath(empty, inpath, 0)
	if out != emptyFileAns {
		t.Errorf("Fail empty input in file mode, want %s, got %s", emptyFileAns, out)
	}
	//check filemode only dir
	fmt.Println("check dir")
	out = checkOutPath(onlydir, inpath, 0)
	if out != onlydirFileAns {
		t.Errorf("Fail only dir input in file mode, want %s, got %s", onlydirFileAns, out)
	}
	//check filemode only filename
	out = checkOutPath(onlyfilename, inpath, 0)
	if out != onlyfilenameFileAns {
		t.Errorf("Fail only filename input in file mode, want %s, got %s", onlyfilenameFileAns, out)
	}
	//check filemode with fullpath
	out = checkOutPath(fullpath, inpath, 0)
	if out != fullpathFileAns {
		t.Errorf("Fail fullpath input in file mode, want %s, got %s", fullpath, out)
	}
	//check dirmode empty
	out = checkOutPath(empty, inpath, 1)
	if out != emptyDirAns {
		t.Errorf("Fail empty input in Dir mode, want %s, got %s", emptyDirAns, out)
	}
	//check dirmode only dir
	out = checkOutPath(onlydir, inpath, 1)
	if out != onlydirDirAns {
		t.Errorf("Fail only Dir input in Dir mode, want %s, got %s", onlydirDirAns, out)
	}
	//check dirmode only filename
	out = checkOutPath(onlyfilename, inpath, 1)
	if out != onlyfilenameDirAns {
		t.Errorf("Fail only filename input in Dir mode, want %s, got %s", onlyfilenameDirAns, out)
	}
	//check dirmode with fullpath
	out = checkOutPath(fullpath, inpath, 1)
	if out != fullpathDirAns {
		t.Errorf("Fail fullpath input in Dir mode, want %s, got %s", fullpathDirAns, out)
	}
}

func TestCounter(t *testing.T) {
	numOFCPU = 4
	las, _ := openLasFile("testdata/sample.las")
	gap := float64(0.1)
	DensityMap, _, _ := Calculate(las, las.Header.MinX, las.Header.MinY, las.Header.MaxX, las.Header.MaxY, gap)
	total := 0
	for _, i := range DensityMap {
		total += i
	}
	if total != las.Header.NumberPoints {
		t.Errorf("total point count is not equal want %v, got %v", las.Header.NumberPoints, total)
	}
}
func TestResult(t *testing.T) {
	numOFCPU = 4
	las, _ := openLasFile("testdata/sample.las")
	gap := float64(1)
	result, _, _ := Calculate(las, las.Header.MinX, las.Header.MinY, las.Header.MaxX, las.Header.MaxY, gap)
	baseline := BaseLine(las)

	size := len(baseline)
	if size != len(result) {
		t.Errorf("size is not same want %v, got %v", size, len(result))
	}

	sum1 := 0
	sum2 := 0
	for i := 0; i < size; i++ {
		sum1 += baseline[i]
		sum2 += result[i]
	}
	if sum1 != sum2 {
		t.Errorf("Total point is not same want %v, got %v", sum1, sum2)
	}
	for i := 0; i < size; i++ {
		if baseline[i] != result[i] {
			t.Errorf("Point in %v is not same want %v, got %v", i, baseline[i], result[i])
			break
		}
	}

}

func BenchmarkCounter(b *testing.B) {
	numOFCPU = 4
	las, _ := openLasFile("testdata/sample.las")
	gap := float64(0.1)
	b.Run("BaseLine", func(b *testing.B) {
		b.ResetTimer()
		BaseLine(las)
	})
	b.Run("Solution1", func(b *testing.B) {
		b.ResetTimer()
		Calculate(las, las.Header.MinX, las.Header.MinY, las.Header.MaxX, las.Header.MaxY, gap)
	})

}

//BaseLine shows the baseline result
func BaseLine(lf *lidario.LasFile) []int {
	c := &counter.Counter{}
	c.Init(lf.Header.MinX, lf.Header.MinY, lf.Header.MaxX, lf.Header.MaxY, 1)
	read := counter.Reader{}
	go read.Serve(c.ReadStream, c)
	write := counter.Writer{}
	go write.Serve(c.WriteStream, c)
	//counter start counting
	wg := &sync.WaitGroup{}
	wg.Add(1)
	go c.Count(0, lf.Header.NumberPoints, lf, wg)
	wg.Wait()
	fmt.Printf("Baseline Total Point %d \n", c.NumberOfPoint)
	return c.DensityMap
}
