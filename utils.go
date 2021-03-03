// Copyright 2019 Hong-Ping Lo. All rights reserved.
// Use of this source code is governed by a BDS-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strings"

	"github.com/jblindsay/lidario"
)

func openLasFile(fileName string) (*lidario.LasFile, error) {
	return lidario.NewLasFile(fileName, "r")
}

func openLasHeader(fileName string) (*lidario.LasFile, error) {
	return lidario.NewLasFile(fileName, "rh")
}

func findFile(root string, match string) (file []string) {
	fmt.Println("Finding Las File in :")
	fmt.Println(root)
	err := filepath.Walk(root, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			fmt.Printf("prevent panic by handling failure accessing a path %q: %v\n", path, err)
			return err
		}
		if info.IsDir() {
			//fmt.Printf("skipping a dir without errors: %+v \n", info.Name())
			return nil
		}

		if strings.Contains(strings.ToLower(info.Name()), match) {
			file = append(file, path)
			return nil
		}
		return nil
	})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("Found", len(file), "las file")
	return file
}
