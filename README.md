# C64Gfx
Utility to convert image files to C64 binary files

Find the current binary in [releases](https://github.com/Sakrac/C64Gfx/releases)

## Usage

C64Gfx [-palette=<image>] -{type} <source image> [additional type params]

where {type} is one of:
 * -spiral: generate ordered spiral fade table
 * -fade2: generate 8-step palette fade table
 * -fadecols: generate 16-step per-color fade table
 * -charspr: split image into char/sprite data
 * -columns: export custom columns
 * -bobfont: export bob font .bin + .wid
 * -agnus: export 2x2 packed grayscale data
 * -textmc: multicolor text picture
 * -bitmap: hires bitmap
 * -bitmapmc: multicolor bitmap
 * -multisprite: export a large sprite cut up into hardware sprites
 * -screens: merge multiple screen/charset sets into one charset
 * -bundle: combine font data for multiple screens into a single font
 * -texthires: hires text picture
 * -rows: convert to 8 pixels per byte row by row
 * -ebcm: extended background color mode (text)
 * \<no arg>: convert EBCM image from 4 background colors


## -columns

C64Gfx -columns \<image> \<out> \<bg> count dim [-mc=col01,col10,col11] [-oc=col] [-oci=col] [-pad=x]

Example: `C64Gfx -columns image.png output 0 5 3x21 -mc=1,2,3 -pad=1`

* image - image file
* out - output file
* bg - background color
* count - number of columns, can be either in a row or a grid by specifying XxY (for example 5x2)
* dim - dimension of each column in [bytes]x[lines] for example 3x21 for a sprite or 1x8 for chars
* -mc=c01,c10,c11 - export multi color binary with listed colors
* -oc=x - write outline mask to a separate output file
* -oci=x - interleave outline mask in output
* -pad - insert x bytes between each column (-pad=1 for sprites)

## -textmc

Text mode multicolor image conversion

C64Gfx -textmc \<image> \<bg> \<col0> \<col1> -out=\<out> [-wid=char width] [-hgt=char height] [-skip0] [-rawcol]

Example: `C64Gfx -textmc picture.png 6 14 3 -out=picture`

* image - image file
* bg - background color ($d020)
* col0 - multicolor 0
* col1 - multicolor 1
* -out=file - output file w/o extension, adds .chr for char data, .scr for screen data, .col for color data
* -skip0 - skip the first character (this is always empty)
* -rawcol - each color in the .col file is 1 byte, by default each byte has 2 colors

## -texthires

C64Gfx -texthires \<image> [-bg=col] [-out=file] [-skip0] [-rawcol]

Example: `C64Gfx -texthires picture.png -out=picture`

Creates .chr, .scr and .col files.

## -bitmap

C64Gfx -bitmap \<image> [-out=\<out>] [-png=\<png>] [-dither=\<1-64>]

Example: `C64Gfx -bitmap picture.png -out=picture -png=picture_preview`

## -bitmapmc

C64Gfx -bitmapmc \<image> \<bg> [-out=\<out>] [-koala=\<koala-file>] [-png=\<png-file>] [-wid=char width] [-hgt=char height] [-rawcol] [-count=num] [-dither=\<1-64>] [-subst=\<subst.png>,col0,col1..]

Example: `C64Gfx -bitmapmc picture.png 0 -koala=picture`

## -ebcm

C64Gfx -ebcm \<image> \<bg0> \<bg1> \<bg2> \<bg3> [-out=\<out>] [-rawcol] [-skip0] [-uniq=\<png>] [-wid=char width] [-hgt=char height]

Example: `C64Gfx -ebcm picture.png 0 1 2 3 -out=picture`

## -screens

C64Gfx -screens \<out.chr> \<in1> [in2] [in3] ...

Example: `C64Gfx -screens output.chr screen1 screen2 screen3`

* uses \<in>.chr + \<in>.scr
* prefix input with * for EBCM input sets

## -bundle

C64Gfx -bundle \<out.chr> \<in1> [in2] [in3] ...

Example: `C64Gfx -bundle output.chr screen1 screen2 screen3`

Reads \<in>.chr and writes \<in>.chrmap files.

## -charspr

C64Gfx -charspr \<image> \<out> \<bg color> \<sprite color>

Example: `C64Gfx -charspr picture.png output 0 1`

## -rows

C64Gfx -rows \<image> \<out> \<col>

Example: `C64Gfx -rows picture.png output 1`

## Default mode (no type switch)

C64Gfx \<image> \<bg0> \<bg1> \<bg2> \<bg3> [-out=\<out>] [-skip0] [-rawcol]

Example: `C64Gfx picture.png 0 1 2 3 -out=picture`

## -multisprite

C64Gfx -multisprite \<png> \<out> \<bg> \<mc1> \<mc0> [-overlay=color] [-singlecount]

Example: `C64Gfx -multisprite sprites.png sprites_out 0 1 2 -overlay=3 -singlecount`

* png - any image file
* out - path & file to output
* bg - color index for background color
* mc1 - multicolor 1 color index
* mc0 - multicolor 0 color index
* -overlay=x - x is color for hi-res overlay sprites

## Build

C64Gfx should compile with most C compilers on most modern platforms without modification.

This repository ships with CMake presets for Clang on Linux and Windows, and for generating Visual Studio projects on Windows.

### Linux (Clang)

```bash
cmake --preset linux-debug
cmake --build --preset linux-debug
```

Release build:

```bash
cmake --preset linux-release
cmake --build --preset linux-release
```

### Windows (Clang + Ninja)

```powershell
cmake --preset windows-clang-debug
cmake --build --preset windows-clang-debug
```

### Windows (Visual Studio project generation)

MSVC toolchain:

```powershell
cmake --preset windows-vs-msvc
cmake --build --preset windows-vs-msvc-debug
```

ClangCL toolchain in Visual Studio generator:

```powershell
cmake --preset windows-vs-clangcl
cmake --build --preset windows-vs-clangcl-debug
```

x86/Win32 (MSVC):

```powershell
cmake --preset windows-vs-msvc-x86
cmake --build --preset windows-vs-msvc-x86-debug
```

x86/Win32 (ClangCL):

```powershell
cmake --preset windows-vs-clangcl-x86
cmake --build --preset windows-vs-clangcl-x86-debug
```

