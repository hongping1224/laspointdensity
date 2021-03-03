// Copyright 2019 Hong-Ping Lo. All rights reserved.
// Use of this source code is governed by a BDS-style
// license that can be found in the LICENSE file.

package counter

import (
	"math"
	"sync"

	"github.com/jblindsay/lidario"
)

type Point struct {
	X float64
	Y float64
}

const (
	allReturn   = 0
	firstReturn = 1
	lastReturn  = 2
)

type Counter struct {
	LBCoordinate    Point
	NumberOfPoint   int
	DensityMap      []int
	MapSize         int
	XRange          int
	YRange          int
	gap             float64
	minZ            float64
	maxZ            float64
	ReadStream      chan Point
	WriteStream     chan int
	WG              *sync.WaitGroup
	CalculationMode int
}

func (counter *Counter) Count(start int, end int, lf *lidario.LasFile, wg *sync.WaitGroup) {
	counter.WG = wg

	switch counter.CalculationMode {

	case firstReturn:
		for i := start; i < end; i++ {
			lp, _ := lf.LasPoint(i)
			if lp.IsFirstReturn() {
				point := lp.PointData()
				if point.Z < counter.minZ || point.Z > counter.maxZ {
					continue
				}
				p := Point{X: point.X, Y: point.Y}
				counter.ReadStream <- p
			}
		}
	case lastReturn:
		for i := start; i < end; i++ {
			lp, _ := lf.LasPoint(i)
			if lp.IsLateReturn() {
				point := lp.PointData()
				if point.Z < counter.minZ || point.Z > counter.maxZ {
					continue
				}
				p := Point{X: point.X, Y: point.Y}
				counter.ReadStream <- p
			}
		}
	default:
		for i := start; i < end; i++ {
			x, y, z, _ := lf.GetXYZ(i)
			if z < counter.minZ || z > counter.maxZ {
				continue
			}
			p := Point{X: x, Y: y}
			counter.ReadStream <- p
		}
	}

	close(counter.ReadStream)
}

func (counter *Counter) Init(minx, miny, maxx, maxy, gap float64, mode int, minZ, maxZ float64) {
	counter.ReadStream = make(chan Point)
	counter.WriteStream = make(chan int)
	//DensityMap

	//LBCoordinate
	counter.LBCoordinate = Point{X: minx, Y: miny}
	counter.gap = gap
	counter.XRange = int(math.Ceil((maxx - minx) / gap))
	counter.YRange = int(math.Ceil((maxy - miny) / gap))
	counter.minZ = minZ
	counter.maxZ = maxZ
	counter.DensityMap = make([]int, counter.XRange*counter.YRange)
	counter.MapSize = counter.XRange * counter.YRange
	counter.CalculationMode = mode

}

type Reader struct {
	done bool
}

func (reader *Reader) Serve(points <-chan Point, counter *Counter) {
	reader.done = false
	for {
		p, more := <-points
		if more {
			readerHandler(p, counter)
		} else {
			close(counter.WriteStream)
			break
		}
	}
	reader.done = true
}

func readerHandler(p Point, counter *Counter) {
	//calculate coordinate of point and  put it into stream
	counter.WriteStream <- xyToIndex(p.X, p.Y, counter)
}

func xyToIndex(x float64, y float64, counter *Counter) int {
	xr := (x - counter.LBCoordinate.X) / counter.gap
	yr := (y - counter.LBCoordinate.Y) / counter.gap
	//fmt.Println(counter.LBCoordinate.X, counter.LBCoordinate.Y)
	index := int(math.Floor(xr)) + (int(math.Floor(yr)) * counter.XRange)
	return index
}

type Writer struct {
}

func (writer *Writer) Serve(indexs <-chan int, counter *Counter) {

	for {
		index, more := <-indexs
		if more {
			writerHandler(index, counter)
		} else {
			break
		}
	}
	counter.WG.Done()

}

func writerHandler(index int, counter *Counter) {
	if index < 0 || index >= counter.MapSize {
		//ignore point outside map
		return
	}
	counter.DensityMap[index]++
	counter.NumberOfPoint++
}
