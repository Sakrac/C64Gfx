# C64Gfx
 Utility to convert image files to C64 binary files

## Build

C64Gfx should compile with most C compilers on most modern platforms without modification.

## Usage

C64Gfx [-palette=<image>] -{type} <source image> [additional type params]

where {type} is one of:
 * -columns: export custom columns (enter without params for info)
 * -textmc: multicolor text picture (enter without params for info)
 * -multisprite: export a large sprite cut up into hardware sprites (enter without params for info)
 * -bundle: combine font data for multiple screens into a single font
 * -texthires: hires text picture (enter without params for info)
 * \<no arg>: convert extended color background image (enter without params for info));


## -columns

C64Gfx -columns \<image> \<out> \<bg> count dim [-mc=col01,col10,col11] [-pad=x]

* image - image file
* out - output file
* bg - background color
* count - number of columns, can be either in a row or a grid by specifying XxY (for example 5x2)
* dim - dimension of each column in [bytes]x[lines] for example 3x21 for a sprite or 1x8 for chars
* -mc=c01,c10,c11 - export multi color binary with listed colors
* -pad - insert x bytes between each column (-pad=1 for sprites)

## -textmc

Text mode multicolor image conversion

C64Gfx -textmc \<image> \<bg> \<col0> \<col1> -out=\<out>

* image - image file
* bg - background color ($d020)
* col0 - multicolor 0
* col1 - multicolor 1
* -out=file - output file w/o extension, adds .chr for char data, .scr for screen data, .col for color data
* -skip0 - skip the first character (this is always empty)
* -rawcol - each color in the .col file is 1 byte, by default each byte has 2 colors

## -texthires

Similar arguments as -textmc

## -screens

This combines multiple exported screens exported from pixcen into a single character set

## -bundle

C64Gfx -bundle \<out> \<in set 1> \<in set 2> \<in set 3>

Combines multiple exported screen differently than -screens

## -charspr

C64Gfx -charspr \<image> \<out> \<bg color> \<sprite color>

## -multisprite

C64Gfx -multisprite \<png> \<out> \<bg> \<mc1> \<mc0> -overlay=color

* png - any image file
* out - path & file to output
* bg - color index for background color
* mc1 - multicolor 1 color index
* mc0 - multicolor 0 color index
* -overlay=x - x is color for hi-res overlay sprites

