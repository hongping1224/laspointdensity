// Copyright 2019 Hong-Ping Lo. All rights reserved.
// Use of this source code is governed by a BDS-style
// license that can be found in the LICENSE file.

package main

import (
	"encoding/binary"
	"flag"
	"fmt"
	"image"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"time"
	"unsafe"

	"github.com/hongping1224/laspointdensity/counter"
	tiff32 "github.com/hongping1224/laspointdensity/tiff"
	"github.com/jblindsay/lidario"
)

var numOFCPU int

func main() {
	numOFCPU = runtime.NumCPU()
	flag.IntVar(&numOFCPU, "cpuCount", numOFCPU, "Cpu use for compute")
	Mode := 0
	flag.IntVar(&Mode, "mode", Mode, "Calculation Mode (default 0): \n 0 for calculate pointcloud seperately \n 1 for calculate combine density")
	file := ""
	flag.StringVar(&file, "file", file, "specific filename to process")
	dir := "./"
	flag.StringVar(&dir, "dir", dir, "directory to process")
	gap := float64(1.0)
	flag.Float64Var(&gap, "size", gap, "pixel size")
	minZ := -1.0
	flag.Float64Var(&minZ, "minZ", minZ, "minZ")
	maxZ := math.MaxFloat64
	flag.Float64Var(&maxZ, "maxZ", maxZ, "maxZ")
	outPath := ""
	flag.StringVar(&outPath, "out", outPath, "Output path (default PointDensity_xxx.tiff)\n single file or Mode 1, Input a filepath \n Mode 0 , Input a Dir")
	returnOption := 0
	flag.IntVar(&returnOption, "return", returnOption, "Point Return Mode (default 0): \n 0 for all return\n 1 for only first return \n 2 for only last return")
	flag.Parse()
	start := time.Now()
	defer exit(start)
	fmt.Printf("Running Program on %v Thread(s)\n", numOFCPU)
	runtime.GOMAXPROCS(numOFCPU)

	if file != "" {
		fmt.Println("Calculating", file)
		outPath = checkOutPath(outPath, file, 0)
		Run(file, gap, outPath, returnOption, minZ, maxZ)
		return
	}
	//run on a directory
	//check directory exist

	if _, err := os.Stat(dir); os.IsNotExist(err) {
		log.Fatal(err)
		return
	}

	//find all las file
	fmt.Println(dir)
	lasfile := findFile(dir, ".las")
	//if mode 0, output each time
	if Mode == 0 {
		for _, path := range lasfile {
			o := checkOutPath(outPath, path, 0)
			fmt.Println("Calculating", path)
			Run(path, gap, o, returnOption, minZ, maxZ)
		}
		return
	}
	//if mode 1, combine output
	minx, miny, maxx, maxy := findBoundary(lasfile)
	totalPoint := totalPoint(lasfile)
	outPath = checkOutPath(outPath, dir, 1)
	mapa := BatchRun(lasfile, minx, miny, maxx, maxy, gap, minZ, maxZ, outPath, returnOption)
	sum := 0
	for _, c := range mapa {
		sum += c
	}
	fmt.Println("total", totalPoint, "sum", sum)
	return

}

func exit(start time.Time) {
	end := time.Now()
	elp := end.Sub(start)
	fmt.Println("Finish Job, Used :", elp)
}

//SaveTiff Save DensityMap into tiff
func SaveTiff(filename string, DensityMap []int, w, h int, gap float64, LB counter.Point) {
	//convert to 2d array
	//index := int(math.Floor(xr)) + (int(math.Floor(yr)) * w)
	fmt.Println("Creating Image")
	im := tiff32.NewGray32(image.Rectangle{Min: image.Point{X: 0, Y: 0}, Max: image.Point{X: w, Y: h}})
	for y := h - 1; y > -1; y-- {
		for x := 0; x < w; x++ {
			i := x + (y * w)
			im.SetGray32(x, h-1-y, tiff32.Gray32Color{Y: uint32(DensityMap[i])})
		}
	}
	//write tfw file
	fmt.Println("Saving Tiff file at:", filename)
	f, err := os.Create(filename + ".tiff")
	if err != nil {
		log.Fatal(err)
		return
	}
	err = tiff32.Encode(f, im, nil)
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
	tfw, err := os.Create(filename + ".tfw")
	gaps := fmt.Sprintf("%g", gap)
	tfw.WriteString(gaps + "\n0\n0\n-" + gaps + "\n")
	tfw.WriteString(fmt.Sprintf("%.5f", LB.X+(gap/2)) + "\n")
	tfw.WriteString(fmt.Sprintf("%.5f", LB.Y+(gap*(float64(h)-0.5))))
	return
}

//SavePointDensityTiff Save DensityMap into tiff
func SavePointDensityTiff(filename string, DensityMap []int, w, h int, gap float64, LB counter.Point) {
	//convert to 2d array
	//index := int(math.Floor(xr)) + (int(math.Floor(yr)) * w)
	fmt.Println("Creating Image")
	area := float32(gap) * float32(gap)
	im := tiff32.NewGrayFloat32(image.Rectangle{Min: image.Point{X: 0, Y: 0}, Max: image.Point{X: w, Y: h}})
	for y := h - 1; y > -1; y-- {
		for x := 0; x < w; x++ {
			i := x + (y * w)

			a := float32(DensityMap[i]) / area
			b := (*[4]byte)(unsafe.Pointer(&a))[:]
			bits := binary.LittleEndian.Uint32(b)
			im.SetGray32(x, h-1-y, tiff32.GrayFloat32Color{Y: bits})

		}
	}
	//write tfw file
	fmt.Println("Saving Tiff file at:", filename)
	f, err := os.Create(filename + "_normalize.tiff")
	if err != nil {
		log.Fatal(err)
		return
	}
	err = tiff32.Encode(f, im, nil)
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
	tfw, err := os.Create(filename + "_normalize.tfw")
	gaps := fmt.Sprintf("%g", gap)
	tfw.WriteString(gaps + "\n0\n0\n-" + gaps + "\n")
	tfw.WriteString(fmt.Sprintf("%.5f", LB.X+(gap/2)) + "\n")
	tfw.WriteString(fmt.Sprintf("%.5f", LB.Y+(gap*(float64(h)-0.5))))
	return
}

//Run Calculate Point Density On a Single File
func Run(filepath string, gap float64, outPath string, returnOption int, minZ, maxZ float64) {

	lf, err := openLasFile(filepath)
	if err != nil {
		log.Println(err)
		return
	}
	fmt.Println(lf.Header)
	densityMap, xrange, yrange := Calculate(lf, lf.Header.MinX, lf.Header.MinY, lf.Header.MaxX, lf.Header.MaxY, gap, returnOption, minZ, maxZ)
	LB := counter.Point{X: lf.Header.MinX, Y: lf.Header.MinY}
	SaveTiff(outPath, densityMap, xrange, yrange, gap, LB)
	SavePointDensityTiff(outPath, densityMap, xrange, yrange, gap, LB)
}

//BatchRun Calculate and Combine Point Density On Many File
func BatchRun(files []string, minx, miny, maxx, maxy, gap, minZ, maxZ float64, outPath string, returnOption int) []int {
	//finalMap
	var densityMap []int
	//w
	var xrange int
	//h
	var yrange int
	for i, path := range files {
		fmt.Println("Calculating(", i, "/", len(files), ")", path)
		las, err := openLasFile(path)

		if err != nil {
			log.Fatal(err)
		}

		var density []int
		density, xrange, yrange = Calculate(las, minx, miny, maxx, maxy, gap, returnOption, minZ, maxZ)
		if densityMap == nil {
			densityMap = density
		} else {
			for i, count := range density {
				densityMap[i] += count
			}
		}
		las.Close()
		runtime.GC()
	}
	LB := counter.Point{X: minx, Y: miny}
	SaveTiff(outPath, densityMap, xrange, yrange, gap, LB)
	SavePointDensityTiff(outPath, densityMap, xrange, yrange, gap, LB)
	return densityMap
}

func totalPoint(files []string) (total int) {
	total = 0
	for _, path := range files {
		las, err := openLasHeader(path)
		if err != nil {
			log.Fatal(err)
		}
		total += las.Header.NumberPoints
		las.Close()
	}
	return
}

func findBoundary(files []string) (minx, miny, maxx, maxy float64) {
	minx = math.MaxFloat64
	miny = math.MaxFloat64
	maxx = 0.0
	maxy = 0.0
	for _, path := range files {
		las, err := openLasHeader(path)
		if err != nil {
			log.Fatal(err)
		}
		minx = math.Min(las.Header.MinX, minx)
		miny = math.Min(las.Header.MinY, miny)
		maxx = math.Max(las.Header.MaxX, maxx)
		maxy = math.Max(las.Header.MaxY, maxy)
		las.Close()
	}
	return
}

//Calculate setup 1 reader and 1 writer for each Counter
//each Counter reads write 1/NumOfCPU of points, and add em up together in the end
func Calculate(lf *lidario.LasFile, MinX, MinY, MaxX, MaxY, gap float64, returnOption int, minZ, maxZ float64) ([]int, int, int) {
	counters := make([]*counter.Counter, numOFCPU)
	//make counter accoding to num CPU
	for i := 0; i < numOFCPU; i++ {
		counters[i] = &counter.Counter{}
		counters[i].Init(MinX, MinY, MaxX, MaxY, gap, returnOption, minZ, maxZ)
		read := counter.Reader{}
		go read.Serve(counters[i].ReadStream, counters[i])
		write := counter.Writer{}
		go write.Serve(counters[i].WriteStream, counters[i])
	}
	//open las file
	//assign stating point, and ending point

	wg := &sync.WaitGroup{}
	wg.Add(numOFCPU)
	nOP := make([]int, numOFCPU+1)
	nOP[0] = 0
	nOP[numOFCPU] = lf.Header.NumberPoints

	for i := 1; i < numOFCPU; i++ {
		nOP[i] = lf.Header.NumberPoints / numOFCPU * i
	}

	for i := 0; i < numOFCPU; i++ {
		go counters[i].Count(nOP[i], nOP[i+1], lf, wg)
	}

	wg.Wait()
	//add em up when all done.
	mapsize := len(counters[0].DensityMap)
	for j := 1; j < numOFCPU; j++ {

		for i := 0; i < mapsize; i++ {
			counters[0].DensityMap[i] += counters[j].DensityMap[i]
		}
	}
	return counters[0].DensityMap, counters[0].XRange, counters[0].YRange
}

//0 for file , 1 for dir
func checkOutPath(outPath, inPath string, mode int) string {
	outdir, outfile := filepath.Split(outPath)
	indir, infile := filepath.Split(inPath)
	fmt.Println(outfile)
	if mode == 0 {
		if outdir == "" {
			outdir = indir
		}
		if outfile == "" {
			outfile = "PointDensity_" + strings.TrimSuffix(infile, ".las")
		} else {
			outfile = strings.TrimSuffix(outfile, filepath.Ext(outfile))
		}
	}
	if mode == 1 {
		outfile = "PointDensity"
		if outdir == "" {
			outdir = indir
		}
	}
	return filepath.Join(outdir, outfile)
}
