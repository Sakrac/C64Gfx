#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#include <stdint.h>
#include <math.h>
#include <stdlib.h>     /* qsort */
#include <windows.h>

// can be overwritten by -palette=<image>
static uint8_t palette[ 16 ][ 3 ] = {
	{ 0, 0, 0 },		// #000000
	{ 255, 255, 255 },	// #FFFFFF
	{ 137, 64, 54 },	// #880000
	{ 122, 191, 199 },	// #AAFFEE
	{ 138, 70, 174 },	// #CC44CC
	{ 104, 169, 65 },	// #00CC55
	{ 62, 49, 162 },	// #0000AA
	{ 208, 220, 113 },	// #EEEE77
	{ 144, 95, 37 },	// #DD8855
	{ 92, 71, 0 },		// #664400
	{ 187, 119, 109 },	// #FF7777
	{ 85, 85, 85 },		// #555555
	{ 128, 128, 128 },	// #808080
	{ 172, 234, 136 },	// #AAFF66
	{ 124, 112, 218 },	// #0088FF
	{ 171, 171, 171 }	// #ABABAB
};

int ColDiff( int colA, int colB )
{
	if( colA == colB ) { return 0; }
	int r = ( int )palette[ colA ][ 0 ] - ( int )palette[ colB ][ 0 ];
	int g = ( int )palette[ colA ][ 1 ] - ( int )palette[ colB ][ 1 ];
	int b = ( int )palette[ colA ][ 2 ] - ( int )palette[ colB ][ 2 ];

	return r*r+g*g+b*b;
}

uint8_t* LoadPicture( const char* file, int* wr, int* hr )
{
	int w, h, n;
	uint8_t *raw = stbi_load( file, &w, &h, &n, 0 ), *src = raw;
	if (!raw) { return 0; }
	uint8_t *img = (uint8_t*)malloc( w * h ), *dst = img;
	if (img) {
		for (int p = 0, np = w * h; p < np; ++p) {
			int r = src[0], g = src[1], b = src[2];
			int bst = (1 << 30) - 1;
			int bc = 0;
			for (int c = 0; c < 16; ++c) {
				int or = r - palette[c][0];
				int og = g - palette[c][1];
				int ob = b - palette[c][2];
				int o = or *or +og * og + ob * ob;
				if (o < bst) {
					bst = o;
					bc = c;
				}
			}
			*dst++ = (uint8_t)bc;
			src += n;
		}
	}
	if( wr ) { *wr = w; }
	if( hr ) { *hr = h; }
	free( raw );
	return img;
}

// R' = R/255
// G' = G/255
// B' = B/255
// 
// Cmax = max(R', G', B')
// Cmin = min(R', G', B')
// 
// delta = Cmax - Cmin
// 
// Hue calculation:
// Max = R => H = mod((G-B)/delta,6.0f)
// Max = G => H = mod((B-R)/delta + 2,6.0f)
// Max = B => H = mod((R-G)/delta + 4,6.0f)
//
// Saturation calculation:
// delta == 0 => S=0 else S = delta / (1-|2*L-1 )
// Lightness calculation:
// L = (Cmax-Cmin)/2
// 
// L = (Cmax + Cmin) / 2

struct float3 { float x, y, z; };
struct float3 RGB2HSL( struct float3 rgb )
{
	float r = rgb.x, g = rgb.y, b = rgb.z;
	float cmax = r > g ? (r > b ? r : g > b ? g : b) : g > b ? g : b;
	float cmin = r < g ? (r < b ? r : g < b ? g : b) : g < b ? g : b;
	float delta = cmax - cmin;
	float lightness = 0.5f * (cmin + cmax);
	float hue = 0.0f;
	float saturation = 0.0f;
	if( delta > 0.0001f ) {
		saturation = delta / (1 - fabsf( 2 * lightness - 1.0f ));
		if( r > g && r > b ) { hue = (g - b) / delta; }
		else if( g > b ) { hue = (b - r) / delta + 2.0f; }
		else { hue = (r - g) / delta + 4.0f; }
		if( hue < 0.0f ) { hue += 6.0f; }
	}
	struct float3 hsl = { hue, saturation, lightness };
	return hsl;
}
struct float3 HSL2RGB( struct float3 hsl )
{
	float c = ( 1 - fabsf( 2 * hsl.z - 1 ) ) * hsl.y;
	float x = c * (1 - fabsf( fmodf( hsl.x, 2.0f ) - 1.0f ) );
	float m = hsl.z - 0.5f*c;
	struct float3 rgb;
	float H = fmodf( hsl.z + 6.0f, 6.0f );
	if( hsl.x < 1 ) { rgb.x = c; rgb.y = x; rgb.z = 0; }
	else if( hsl.x < 2 ) { rgb.x = x; rgb.y = c; rgb.z = 0; }
	else if( hsl.x < 3 ) { rgb.x = 0; rgb.y = c; rgb.z = x; }
	else if( hsl.x < 4 ) { rgb.x = 0; rgb.y = x; rgb.z = c; }
	else if( hsl.x < 5 ) { rgb.x = x; rgb.y = 0; rgb.z = c; }
	else { rgb.x = c; rgb.y = 0; rgb.z = x; }
	struct float3 ret = { rgb.x + m, rgb.y + m, rgb.z + m };
	return ret;
}

float HSLDist( struct float3 hsl1, struct float3 hsl2 )
{

	float h1 = hsl1.x, s1 = hsl1.y, l1 = hsl1.z, h2 = hsl2.x, s2 = hsl2.y, l2 = hsl2.z;

	float hd = h2 - h1;
	if( hd > 3.0f ) { hd -= 6.0f; }
	else if( hd < -3.0f ) { hd += 6.0f; }
	float sd = s2 - s1;
	float ld = l2 - l1;

	return sqrtf( hd*hd / 9.0f + sd*sd + ld*ld );
}

uint8_t Dither8x8[8][9] = {
	{ 0, 48, 12, 60, 3, 51, 15, 63 },
	{ 32, 16, 44, 28, 35, 19, 47, 31 },
	{ 8, 56, 4, 52, 11, 59, 7, 55 },
	{ 40, 24, 36, 20, 43, 27, 39, 23 },
	{ 2, 50, 14, 62, 1, 49, 13, 61 },
	{ 34, 18, 46, 30, 33, 17, 45, 29 },
	{ 10, 58, 6, 54, 9, 57, 5, 53 },
	{ 42, 26, 38, 22, 41, 25, 37, 21 } };


char* GetSwitch( const char* match, char** swtc, int swtn )
{
	int l = ( int )strlen( match );
	while( swtn )
	{
		if( _strnicmp( match, *swtc, l ) == 0 )
		{
			if( ( *swtc )[ l ] == '=' ) return *swtc + l + 1;
			else return *swtc;
		}
		++swtc;
		--swtn;
	}
	return 0;
}

float square( float x ) { return x * x; }

#define MAX_ARGS 32
#define INV_255 ( 1.0f / 255.0f )

#define NUM_SPIRAL_POINTS 512
#define M_PIF       3.14159265358979323846f

struct DistIdx {
	float dist;
	int idx;
};

int sortDistIdx( const void * a, const void * b )
{
	struct DistIdx a32 = *(struct DistIdx* )a;
	struct DistIdx b32 = *(struct DistIdx* )b;
	if( a32.dist == b32.dist ) { return 0; }
	else if( a32.dist < b32.dist ) { return -1; }
	return 1;
}


const int sprw = 24;
const int sprh = 21;
int FindSprite(uint8_t* img, int w, int h, uint8_t bg, uint8_t col, uint8_t mc1, uint8_t mc2, uint8_t* out, int* offs )
{
	uint8_t* pix = img;

	for (int y = 0; y<h; ++y) {
		for (int x = 0; x<w; ++x) {
			uint8_t c = img[x+y*w];
			if (c==col||c==mc1||c==mc2) {
				// find left/right of this color in image
				uint8_t* row = img + y * w;
				int l = x, ls=x;
				int r = x, rs=x;
				int sh = (y<(h-sprh)) ? sprh : (h-y);
				for (int ys = 0; ys<sh; ++ys) {
					for (int xs = 0; xs<w; ++xs) {
						uint8_t c2 = *row++;
						if (c2==col||c2==mc1||c2==mc2) {
							l = l<xs ? l : xs;
							r = r>xs ? r : xs;
							if (xs<ls&&(rs-xs)<=sprw) {
								ls = xs; }
							else if (xs>rs&&(xs-ls)<=sprw) { rs = xs; }
						}
					}
				}
				// for now do ls->rs until better algo..
				if (mc1<16) { // multicolor?
					ls = l;
					rs = (ls+24<w) ? (ls+23) : (w-1);
					if ((rs-ls+1)>sprw) { rs = ls+sprw-1; }
					ls &= 0xfe;
					if (ls>(w-24)) { ls = w-24; }
					rs = (rs+1)&0xfe; // this might make it 26 wide..
					if ((rs-ls)>24) { rs = ls+24; }
					if (col>=16) { // check if there is a sprite color
						for (int ys = 0; ys<sh&&col>=16; ++ys) {
							uint8_t* spr = img+ls+(y+ys) * w;
							for (int xs = 0; xs<(rs-ls); ++xs) {
								uint8_t c2 = *spr++;
								if (c2!=bg && c2!=mc1 && c2!=mc2) {
									col = c2;
									break;
								}
							}
						}
					}
					{
						offs[0] = ls;
						offs[1] = y;
						offs[2] = col;
						int sw = rs-ls;
						for (int ys = 0; ys<sh; ++ys) {
							uint8_t b = 0;
							for (int xs = 0; xs<sw; xs += 2) {
								uint8_t* spr = img+ls+xs+(y+ys) * w;
								uint8_t c2 = *spr;
								if (c2==bg) { c2 = spr[1]; }
								b <<= 2;
								if (c2==col) { b |= 2; } else if (c2==mc1) { b |= 1; } else if (c2==mc2) { b |= 3; }
								if (c2==col||c2==mc1||c2==mc2) {
									spr[0] = bg; spr[1]=bg;
								}
								if ((xs&6)==6) { *out++ = b; }
							}
						}
						return 1;
					}
				} else {
					ls = l;
					rs = (ls+24<w) ? (ls+23) : (w-1);
					if ((rs-ls+1)>sprw) { rs = ls+sprw-1; }
					offs[0] = ls;
					offs[1] = y;
					offs[2] = col;
					int sw = rs+1-ls;
					for (int ys = 0; ys<sh; ++ys) {
						uint8_t b = 0;
						for (int xs = 0; xs<sw; ++xs) {
							uint8_t* spr = img+ls+xs+(y+ys) * w;
							uint8_t c2 = spr[0];
							b <<= 1;
							if (c2==col) { b |= 1;
								spr[0] = bg; }
							if ((xs&7)==7) { *out++ = b; }
						}
					}
					return 1;
				}
			}
		}
	}
	return 0;
}

#define MAX_MULTISPRITE 32
int MultiSprite(const char** args, int argn, char** swtc, int swtn)
{
	if (argn<6) { printf("Usage:\nC64Gfx -multisprite <png> <out> <bg> <mc1> <mc0> -overlay=color\n"); return 0; }

	int w, h;
	uint8_t* img = LoadPicture(args[1], &w, &h);
	if (!img) {
		printf("Can't open file %s as an image\n", args[1]);
		return 1;
	}
	uint8_t spriteData[MAX_MULTISPRITE][64] = { 0 };
	int spriteOffs[MAX_MULTISPRITE][3] = { 0 }; // 3rd is color
	int nSprites = 0;

	uint8_t cols[3] = { (uint8_t)atoi(args[3]), (uint8_t)atoi(args[4]), (uint8_t)atoi(args[5]) };
	// check for overlay sprites if specified
	const char* ovlstr = GetSwitch("overlay", swtc, swtn);
	if (ovlstr) {
		uint8_t ovlc = atoi(ovlstr);
		printf("Finding sprites in %s (%d, %d) with an overlay color of %d and multicolors %d and %d\n", args[1], w, h, ovlc, cols[1], cols[2]);
		int found = 1;
		while (found) {
			found = 0;
			if (FindSprite(img, w, h, cols[0], ovlc, 0xff, 0xff, spriteData[nSprites], &spriteOffs[nSprites][0])) {
				found = 1;
				nSprites++;
				if (nSprites>=MAX_MULTISPRITE) { break; }
			}
		}
	} else {
		printf("Finding sprites in %s (%d, %d) with multicolors %d and %d\n", args[1], w, h, cols[1], cols[2]);
	}
	int nOutline = nSprites;
	while( nSprites < MAX_MULTISPRITE ) {
		if (FindSprite(img, w, h, cols[0], 0xff, cols[1], cols[2], spriteData[nSprites], &spriteOffs[nSprites][0])) {
			nSprites++;
			if (nSprites>=MAX_MULTISPRITE) { break; }
		} else {
			break;
		}
	}

	char outfile[_MAX_PATH];
	size_t outnameLen = strlen(args[2]);
	memcpy(outfile, args[2], outnameLen);
	strcpy(outfile+outnameLen, ".bin");
	FILE* f = fopen(outfile, "wb");
	if (f) {
		fwrite(spriteData, 64, nSprites, f);
		fclose(f);
	} else {
		printf("could not open sprite output file %s for writing\n", outfile);
		return 1;
	}
	strcpy(outfile+outnameLen, ".s");
	f = fopen(outfile, "w");
	if (f) {
		fprintf(f, "\tdc.b %d\t;single color sprite count\n", nOutline );
		for (int o = 0; o<nOutline; ++o) {
			fprintf(f, "\tdc.b %d, %d, %d\t; offs x, y and sprite color\n", spriteOffs[o][0], spriteOffs[o][1], spriteOffs[o][2]);
		}
		fprintf(f, "\tdc.b %d\t;multicolor sprite count\n", nSprites-nOutline);
		for (int o = nOutline; o<nSprites; ++o) {
			fprintf(f, "\tdc.b %d, %d, %d\t; offs x, y and sprite color\n", spriteOffs[o][0], spriteOffs[o][1], spriteOffs[o][2]);
		}
		// enable sprite mask, multicolor mask
		fprintf(f, "\tdc.b $%02x, $%02x\t; sprite enable mask and multicolor mask\n", (1<<nSprites)-1, ((1<<(nSprites-nOutline))-1)<<nOutline);
		fclose(f);
	} else {
		printf("could not open sprite offset file %s for writing\n", outfile);
		return 1;
	}
	return 0;
}

int LoadPalette(const char* paletteFilename)
{
	int x, y, n;
	uint8_t* data = stbi_load(paletteFilename, &x, &y, &n, 0);
	if (data) {
		uint8_t* scan = data;
		size_t colIdx = 0;
		size_t pixels = x * y, curr = 1;
		for (size_t copy = 0; copy < 3; ++copy) { palette[0][copy] = scan[copy]; }
		while (colIdx < 15 && curr < pixels) {
			size_t comp = 0;
			for (; comp < 3; ++comp) {
				if (palette[colIdx][comp] != scan[comp]) { break; }
			}
			if (comp < 3) {
				++colIdx;
				for (size_t copy = 0; copy < 3; ++copy) { palette[colIdx][copy] = scan[copy]; }
			}
			scan += n;
			curr++;
		}
		stbi_image_free(data);
		return 1;
	}
	return 0;
}

int main( int argc, char* argv[] )
{
	const char* args[MAX_ARGS];
	char* swtc[MAX_ARGS];
	int argn = 0;
	int swtn = 0;
	for( int a = 0; a < argc; ++a ) {
		if( argv[ a ][ 0 ] == '-' ) { swtc[ swtn++ ] = argv[ a ] + 1; }
		else { args[ argn++ ] = argv[ a ]; }
	}

	char currDir[512];

	if (GetSwitch("dir",swtc,swtn)) {
		GetCurrentDirectory(sizeof(currDir), currDir);
		printf("%s\n", currDir);
	}

	if (swtn < 1) {
		printf(
			"gfx [-palette=<image>] -{type} <source image> [additional type params]\n"
			"types:\n"
			" * -spiral: custom demo effect\n"
			" * -fade2: custom demo effect\n"
			" * -fadecols: custom demo effect\n"
			" * -charspr: i don't quite remember\n"
			" * -columns: export custom columns (enter without params for info)\n"
			" * -bobfont: custom demo format font\n"
			" * -agnus: custom demo effect\n"
			" * -textmc: multicolor text picture (enter without params for info)\n"
			" * -multisprite: export a large sprite cut up into hardware sprites (enter without params for info)\n"
			" * -screens: i don't remember\n"
			" * -bundle: combine font data for multiple screens into a single font\n"
			" * -texthires: hires text picture (enter without params for info)\n"
			" * <no arg>: convert extended color background image (enter without params for info)\n");
		return 0;
	}

	const char* paletteFile = GetSwitch("palette", swtc, swtn);
	if (paletteFile) {
		if (!LoadPalette(paletteFile)) {
			printf("Could not open palette image \"%s\"\n", paletteFile);
			return 1;
		}
	}

	if( GetSwitch( "spiral", swtc, swtn ) ) {
		float points[ NUM_SPIRAL_POINTS ][ 2 ];
		float degrees = 720;
		float distLength = 64.0f; // time per pixel away from closest point on curve
		int width = 24;
		int height = 17;
		int length = 100;
		const char* arg;
		arg = GetSwitch( "degrees", swtc, swtn );
		if( arg ) { degrees = (float)atof( arg ); }
		arg = GetSwitch( "distLength", swtc, swtn );
		if( arg ) { distLength = ( float )atof( arg ); }
		arg = GetSwitch( "width", swtc, swtn );
		if( arg ) { width = atoi( arg ); }
		arg = GetSwitch( "height", swtc, swtn );
		if( arg ) { height = atoi( arg ); }
		arg = GetSwitch( "length", swtc, swtn );
		if( arg ) { length = atoi( arg ); }

		float radiusX = (float)(width * 4);
		float radiusY = (float)(height * 4);
		arg = GetSwitch( "radiusX", swtc, swtn );
		if( arg ) { radiusX = ( float )atof( arg ); }
		arg = GetSwitch( "radiusY", swtc, swtn );
		if( arg ) { radiusY = ( float )atof( arg ); }
		for( int p=0; p < NUM_SPIRAL_POINTS; ++p ) {
			float t = (float)p / ( float )NUM_SPIRAL_POINTS;
			float a = t * degrees * M_PIF / 180.0f;
			float rx = radiusX * t;
			float ry = radiusY * t;
			points[p][0] = rx * cosf( a ) + radiusX;
			points[p][1] = ry * sinf( a ) + radiusY;
		}
		// sorted by order...
		struct DistIdx* order = ( struct DistIdx* )malloc( width * height * sizeof( struct DistIdx ) );
		// find the closest time distance for each center of the grid
		int num = 0;
		for( int y = 0; y < height; ++y ) {
			float posY = y*8.0f + 4.0f;
			for( int x = 0; x < width; ++x ) {
				float posX = x*8.0f + 4.0f;
				float closestT = 1e36f;

				for( int p=0; p < NUM_SPIRAL_POINTS; ++p ) {
					float dx = posX - points[ p ][ 0 ];
					float dy = posY - points[ p ][ 1 ];
					float l = sqrtf(dx*dx + dy*dy);
					float tSum = p * length / (float)NUM_SPIRAL_POINTS + distLength * l;
					if( tSum < closestT ) {
						closestT = tSum;
					}
				}
				order[ num ].dist = closestT;
				order[ num ].idx = y * 40 + x;
				++num;
			}
		}
		qsort( order, num, sizeof( struct DistIdx ), sortDistIdx );

		printf( "\n; ordered:\n\tdc.w " );
		for( int i = 0, n = width * height; i<n; ++i ) {
			printf( "$%04x + FadeScreen", order[ i ].idx & 0xffff );
			if( (~i)&0x7 ) { printf(", "); }
			else { printf( "\n\tdc.w " ); }
		}

		return 0;
	}

	if( GetSwitch( "fade2", swtc, swtn ) )
	{
		const float graySteps[] = { 0.48f, 0.55f, 0.8f };
		for( int s=0; s<8; ++s ) {
			for( int c = 0; c < 16; ++c ) {
				float r = palette[ c ][0] * INV_255;
				float g = palette[ c ][1] * INV_255;
				float b = palette[ c ][2] * INV_255;
				float v = s<3 ? graySteps[s] : 1.0f;
				float k = s>=3 ? (7-s)/(7.0f-3.0f) : 1.0f;

				float m = ( ( 0.3f * r ) + ( 0.59f * g ) + ( 0.11f * b ) );
				r = k * ( v * m + ( 1 - v ) * r );
				g = k * ( v * m + ( 1 - v ) * g );
				b = k * ( v * m + ( 1 - v ) * b );
				float closest = 1e36f;
				int col = 11;
				for( int cc = 0; cc<16; ++cc ) {

					float rd = r - palette[ cc ][ 0 ] * INV_255;
					float gd = g - palette[ cc ][ 1 ] * INV_255;
					float bd = b - palette[ cc ][ 2 ] * INV_255;
					float d = rd*rd+gd*gd+bd*bd;
					if( d < closest ) {
						closest = d;
						col = cc;
					}
				}
				printf( "$%02x, ", col );
			}
			printf("\n");
		}
		return 0;
	}


	if( GetSwitch( "fadecols", swtc, swtn ) )
	{
		const float inv255 = 1.0f / 255.0f;
		uint8_t fades[ 16 ][ 16 ];
		const float hiL = 0.25f;
		const float hiS = 0.5f;

		for( int c = 0; c < 16; ++c ) {
			struct float3 orig = { inv255 * palette[ c ][ 0 ], inv255 * palette[ c ][ 1 ], inv255 * palette[ c ][ 2 ] };
			struct float3 hsl_orig = RGB2HSL( orig );
			printf( "\tdc.b" );
			for( int f = 0; f < 16; ++f ) {
				// 0 - black, 15 = c
				struct float3 hsl = hsl_orig;
#if 1
				float ff = (f + 1) / 16.0f;
				if( ff < hiS ) { hsl.y = ff / hiS; }
				else {
					float fr = ff - hiS;
					float fs = fr / (1.0f - hiS);
					hsl.y = (1.0f - fs ) + fs * hsl.y;
				}
				//hsl.y *= ff;
				hsl.x = fmodf( hsl.x + 6.0f * (1.0f - ff), 6.0f );
				if( ff < hiL ) { hsl.z = 0.75f * ff / hiL; }
				else { 
					float fr = ff - hiL;
					float fe = 1.0f - fr / (1.0f - hiL);
					hsl.z = 0.75f * fe + ( 1.0f - fe ) * hsl.z;
				}
#else
				if( f < 8 ) {
					hsl = hsl_orig;
					hsl.x = fmodf( hsl_orig.x + 3.0f, 6.0f );
					hsl.y = (1.0f / 8.0f) * f;
					hsl.z = (1.0f / 8.0f) * f;
				} else {
					float t = (f - 8) / 7.0f;
					hsl.x = hsl_orig.x;
					hsl.y = (1.0f - t) + t * hsl_orig.y;
					hsl.z = (1.0f - t) + t * hsl_orig.z;
				}
#endif
				float closest = 1e36f;
				int b = 0;
				struct float3 rgb = HSL2RGB( hsl );
				for( int m = 0; m < 16; ++m ) {
					struct float3 prgb = { inv255 * palette[ m ][ 0 ], inv255 * palette[ m ][ 1 ], inv255 * palette[ m ][ 2 ] };

					float d = square(prgb.x - rgb.x) + square(prgb.y - rgb.y) + square(prgb.z - rgb.z);

					if( d < closest ) {
						closest = d;
						b = m;
					}
				}
				fades[ c ][ f ] = b;
				if( f ) { printf( "," ); }
				printf( " $%02x", b );
			}
			printf( "\n" );
		}
		return 0;
	}

	if( GetSwitch( "charspr", swtc, swtn ) )
	{
		if( argn < 4 ) { printf( "Usage:\nGfx -charspr <image> <out> <bg color> <sprite color>\n" ); return 0; }
		int w, h;
		uint8_t* img = LoadPicture( args[ 1 ], &w, &h );
		if( !img ) { printf("Failed to load image %s\n", args[1] ); return 1; }
		uint8_t bg = ( uint8_t )atoi( args[3] );
		uint8_t fg = ( uint8_t )atoi( args[ 4 ] );
		int sprWide = w / 24;
		int sprHgt = h / 21;
		int chrWide = sprWide * 3;
		int chrHgt = ( sprHgt * 21 + 7 ) / 8;
		uint8_t* sprOut = (uint8_t*)calloc(1, 64 * sprWide * sprHgt );
		uint8_t* chrOut = (uint8_t*)calloc(1, 8 * chrWide * chrHgt );
		uint8_t* chrCol = (uint8_t*)calloc(1,chrWide*chrHgt);

		for( int y = 0; y < chrHgt; ++y ) {
			for( int x = 0; x < chrWide; ++x ) {
				int bestCol = 12;
				int bestDist = 1<<30;
				uint8_t* chr = img + y * 8 * w + x * 8;
				for( int c = 0; c < 16; ++c ) {
					if( c != bg && c != fg ) {
						int dist = 0;
						for( int py = 0; py < 8; ++py ) {
							for( int px = 0; px < 8; ++px ) {
								uint8_t p = ( (x*8+px) < w && (y*8+py) < h ) ? chr[ py * w + px ] : bg;
								if( p != bg && p != fg ) {
									dist += ColDiff( p, c );
								}
							}
						}
						if( dist < bestDist ) {
							bestCol = c;
							bestDist = dist;
						}
					}
				}
				chrCol[x + y * chrWide] = bestCol;
			}
		}
		// set up the char data
		for( int y = 0; y < chrHgt; ++y ) {
			for( int x = 0; x < chrWide; ++x ) {
				uint8_t mid = chrCol[x+y*chrWide];
				uint8_t* chr = img + y * 8 * w + x * 8;
				for( int py = 0; py < 8; ++py ) {
					for( int px = 0; px < 8; ++px ) {
						uint8_t p = ( ( x * 8 + px ) < w && ( y * 8 + py ) < h ) ? chr[ py * w + px ] : bg;
						int distBG = ColDiff(p, bg);
						int distMid = ColDiff(p, mid);
						int distFG = ColDiff(p, fg);
						if( distMid < distBG || distFG < distBG ) {
							chrOut[8 * ( x + y * chrWide ) + py ] |= 1<<(px^7);
						}
					}
				}
			}
		}
		// set up the sprite data
		for( int y = 0; y < sprHgt; ++y ) {
			for( int x = 0; x < sprWide; ++x ) {
				uint8_t* spr = img + y * 21 * w + x * 24;
				for( int py = 0; py < 21; ++py ) {
					for( int px = 0; px < 24; ++px ) {
						uint8_t p = ( ( x * 24 + px ) < w && ( y * 21 + py ) < h ) ? spr[ py * w + px ] : bg;
						int chId = ((y*21+py)>>3) * chrWide + (( x*24+px )>>3);
						uint8_t mid = chrCol[chId];
						int distBG = ColDiff( p, bg );
						int distMid = ColDiff( p, mid );
						int distFG = ColDiff( p, fg );
						if( distFG < distMid && distFG < distBG ) {
							sprOut[ ( x + y * sprWide ) * 64 + py * 3 + (px>>3) ] |= 1<<((px&7)^7);
						}
					}
				}
			}
		}

		size_t outLen = strlen( args[2] );
		char file[ _MAX_PATH ];
		const char* extChr = ".chr";
		const char* extCol = ".col";
		const char* extSpr = ".spr";
		memcpy( file, args[ 2 ], outLen );
		memcpy( file + outLen, extChr, sizeof( extChr ) + 1 );
		FILE* f;
		if( fopen_s( &f, file, "wb" ) == 0 ) {
			fwrite( chrOut, 8*chrWide * chrHgt, 1, f );
			fclose(f);
		}
		memcpy( file + outLen, extCol, sizeof( extCol ) + 1 );
		if( fopen_s( &f, file, "wb" ) == 0 ) {
			fwrite( chrCol, chrWide * chrHgt, 1, f );
			fclose( f );
		}
		memcpy( file + outLen, extSpr, sizeof( extSpr ) + 1 );
		if( fopen_s( &f, file, "wb" ) == 0 ) {
			fwrite( sprOut, sprWide*sprHgt*64, 1, f );
			fclose( f );
		}


		return 0;
	}

	if( GetSwitch( "columns", swtc, swtn ) )
	{
		if( argn < 5 ) { printf( "Usage:\nGfx -columns <image> <out> <bg> count dim [-mc=col01,col10,col11]\n" ); return 0; }
		int w, h;
		uint8_t* img = LoadPicture( args[ 1 ], &w, &h );

		int wc = w / 8;
		int hc = h / 8;

		char *endStr;
		uint8_t bg = ( uint8_t )strtoul( args[ 3 ], &endStr, 10 );
		uint8_t mc[3] = { 1, 2, 3 };


		int pad = 0;
		const char* padStr = GetSwitch( "pad", swtc, swtn );
		if( padStr ) { pad = atoi( padStr ); }

		int countX = 1;
		int countY = 1;
		const char* cntStr = args[ 4 ];
		countX = strtoul( cntStr, &endStr, 10 );
		if( endStr[ 0 ] == 'x' || endStr[ 0 ] == 'X' ) {
			++endStr;
			countY = strtoul( endStr, &endStr, 10 );
		}

		char* mcstr = GetSwitch("mc", swtc, swtn);
		if (mcstr) {
			mc[0] = (uint8_t)strtoul(mcstr, &endStr, 10);
			if (*endStr) {
				mc[1] = (uint8_t)strtoul(++endStr, &endStr, 10);
				if (*endStr) {
					mc[2] = (uint8_t)strtoul(++endStr, &endStr, 10);
				}
			}
		}

		int dimX = 1;
		int dimY = 16;
		const char* dimStr = args[ 5 ];
		dimX = strtoul( dimStr, &endStr, 10 );
		if( endStr[ 0 ] == 'x' || endStr[ 0 ] == 'X' ) {
			++endStr;
			dimY = strtoul( endStr, &endStr, 10 );
		}

		if( wc < ( countX * dimX ) ) { countX = wc / dimX; }
		if( h < ( countY * dimY ) ) { countY = h / dimY; }

		uint8_t* buf = (uint8_t*)malloc(dimX * dimY * countX * countY + pad * countX * countY ), *out = buf;

		for( int row = 0; row < countY; ++row ) {
			for( int col = 0; col < countX; ++col ) {
				uint8_t* cell = img + row * dimY * w + col * dimX * 8;
				for( int y = 0; y < dimY; ++y ) {
					for( int xb = 0; xb < dimX; ++xb ) {
						uint8_t b = 0;
						if (mcstr) {
							for (int x = 0; x<8; x+=2) {
								b <<= 2;
								if (cell[x]==mc[0]) { b |= 1; }
								else if (cell[x]==mc[1]) { b |= 2; }
								else if (cell[x]==mc[2]) { b |= 3; }
							}
						} else {
							for (int x = 0; x<8; ++x) {
								b <<= 1;
								if (cell[x]!=bg) { b |= 1; }
							}
						}
						*out++ = b;
						cell += 8;
					}
					cell += w - dimX*8;
				}
				for( int p = 0; p < pad; ++p ) { *out++ = 0; }
			}
		}
		FILE* f = fopen( args[ 2 ], "wb" );
		if( f ) {
			fwrite( buf, out - buf, 1, f );
			fclose( f );
		} else { printf( "failed to open file \"%s\" for writing\n", args[2] ); return 1; }
		free( buf );
		return 0;
	}

	if( GetSwitch( "bobfont", swtc, swtn ) )
	{
		if( argn < 3 ) { printf( "Usage:\nGfx -bobfont <src> <out> <chars>\n" ); return 0; }
		int w, h;
		uint8_t* img = LoadPicture( args[ 1 ], &w, &h );
		if( !img ) { printf( "could not open %s\n", args[ 1 ] ); return 1; }
		int chars = atoi( args[ 3 ] );
		if( !chars ) { printf( "no chars?\n" ); return 2; }
		uint8_t widths[ 256 ];
		uint8_t buf[ 256 ][ 8 ];
		int wc = w / 8;
		for( int c = 0; c < chars; ++c ) {
			uint8_t* chr = img + (c / wc) * 8 * w + (c %wc) * 8;
			int last = 0;
			for( int y = 0; y < 8; ++y ) {
				uint8_t byte = 0;
				for( int x = 0; x < 5; ++x ) {
					uint8_t bit = chr[ y + x * w ];
					if( bit ) { byte |= 1 << x; }
				}
				if( byte ) { last = y; }
				buf[ c ][ y ] = byte;
			}
			widths[ c ] = last + 1;
		}
		FILE *f;
		char name[ _MAX_PATH ];
		strcpy( strcpy( name, args[ 2 ] ) + strlen( args[ 2 ]), ".bin" );
		if( fopen_s( &f, name, "wb" ) == 0 ) {
			fwrite( buf, 8, chars, f );
			fclose( f );
		} else { printf( "failed to write output to %s\n", name ); return 3; }
		strcpy( strcpy( name, args[ 2 ] ) + strlen( args[ 2 ] ), ".wid" );
		if( fopen_s( &f, name, "wb" ) == 0 ) {
			fwrite( widths, chars, 1, f );
			fclose( f );
		} else { printf( "failed to write output to %s\n", name ); return 4; }
		return 0;
	}

	if( GetSwitch( "agnus", swtc, swtn ) )
	{
		if( argn < 3 ) { printf( "Usage:\nGfx -agnus <image> <out>\n" ); return 0; }
		int w, h, n;
		uint8_t *raw = stbi_load( args[1], &w, &h, &n, 0 ), *src = raw;
		if( raw ) {
			uint8_t* outbuf = (uint8_t*)malloc( (w / 2) * (h / 2) ), *out = outbuf;
			for( int y = 0, ny = h / 2; y < ny; ++y ) {
				for( int x = 0, nx = w / 2; x < nx; ++x ) {
					uint8_t c[ 4 ] = {
						raw[ n * (y * w * 2 + x * 2) ],
						raw[ n * (y * w * 2 + x * 2 + 1) ],
						raw[ n * (y * w * 2 + x * 2 + w) ],
						raw[ n * (y * w * 2 + x * 2 + w + 1) ] };

					for( int ci = 0; ci < 4; ++ci ) {
						if( c[ ci ] > 130 ) { c[ ci ] = 3; }
						else if( c[ ci ] > 68 ) { c[ ci ] = 2; }
						else if( c[ ci ] > 47 ) { c[ ci ] = 1; }
						else { c[ ci ] = 0; }
					}
					*out++ = c[ 0 ] | (c[ 1 ] << 2) | (c[ 2 ] << 4) | (c[ 3 ] << 6);
				}
			}
			free( raw );
			FILE* f = fopen( args[2], "wb" );
			if( f ) {
				fwrite( outbuf, out-outbuf, 1, f );
				fclose( f );
			}
			free( outbuf );
			return 0;
		} else {
			printf( "failed to open file %s\n", args[ 1 ] );
			return 1;
		}

	}

	if( GetSwitch( "textmc", swtc, swtn ) ) {
		if( argn < 5 ) { printf( "Usage:\nGfx -textmc <image> <bg> <col0> <col1> -out=<out>\n" ); return 0; }

		uint8_t cols[ 3 ] = { (uint8_t)atoi( args[ 2 ] ), (uint8_t)atoi( args[ 3 ] ), (uint8_t)atoi( args[ 4 ] ) };

		int w, h;
		uint8_t* img = LoadPicture( args[ 1 ], &w, &h );

		int wc = w / 8;
		int hc = h / 8;

		uint8_t* screen = (uint8_t*)malloc( wc * hc );
		uint8_t* chars = (uint8_t*)malloc( 256 * 8 );
		uint8_t* color = (uint8_t*)malloc( wc * hc );

		int nunChars = 1;	// first char is empty!
		memset( chars, 0, 8 );

		int prevChrCol = 0;

		// prioritize background color 0, then 1 etc.
		for( int y = 0; y < hc; ++y ) {
			for( int x = 0; x < wc; ++x ) {
				uint8_t* s = img + y * 8 * w + x * 8;
				int chrCol = -1;
				for( int yp = 0; yp < 8; ++yp ) {
					for( int xp = 0; xp < 8; xp+=2 ) {
						uint8_t c = s[ xp + yp * w ];
						if( c < 8 && c != cols[ 0 ] && c != cols[ 1 ] && c != cols[ 2 ] ) {
							chrCol = c;
							break;
						}
					}
				}
				color[ y * wc + x ] = chrCol >= 0 ? ( chrCol | 8 ) : ( prevChrCol | 8 );
				if( chrCol > 0 ) { prevChrCol = chrCol; }
				uint8_t chr[ 8 ] = { 0 };
				for( int yp = 0; yp < 8; ++yp ) {
					uint8_t b = 0;
					for( int xp = 0; xp < 8; xp += 2 ) {
						uint8_t c = s[ xp + yp * w ];
						b <<= 2;
						if( c == chrCol ) { b |= 3; }
						else if( c == cols[ 1 ] ) { b |= 1; }
						else if( c == cols[ 2 ] ) { b |= 2; }
					}
					chr[ yp ] = b;
				}
				int c = 0;
				for( ; c < nunChars; ++c ) {
					int same = 1;
					for( int b = 0; b<8; ++b ) {
						if( chars[ c * 8 + b ] != chr[ b ] ) {
							same = 0;
							break;
						}
					}
					if( same ) { break; }
				}
				if( c == nunChars ) {
					if( nunChars < 256 ) {
						for( int b = 0; b<8; ++b ) { chars[ c * 8 + b ] = chr[ b ]; }
						++nunChars;
					}
					else { c = 0; }
				}
				screen[ x + y * wc ] = c;
			}
		}
		printf( "Used chars for %s = %d\n", args[1], nunChars );

		const char* out = GetSwitch( "out", swtc, swtn );
		if( out ) {
			size_t outLen = strlen( out );
			char file[ _MAX_PATH ];
			const char* extChr = ".chr";
			const char* extScr = ".scr";
			const char* extCol = ".col";
			memcpy( file, out, outLen );
			memcpy( file + outLen, extChr, sizeof( extChr ) + 1 );
			FILE* f = fopen( file, "wb" );
			if( f ) {
				if( GetSwitch( "skip0", swtc, swtn ) ) {
					fwrite( chars + 8, nunChars * 8 - 8, 1, f );
				}
				else {
					fwrite( chars, nunChars * 8, 1, f );
				}
				fclose( f );
			}

			memcpy( file, out, outLen );
			memcpy( file + outLen, extScr, sizeof( extScr ) + 1 );
			f = fopen( file, "wb" );
			if( f ) {
				fwrite( screen, wc*hc, 1, f );
				fclose( f );
			}

			memcpy( file, out, outLen );
			memcpy( file + outLen, extCol, sizeof( extCol ) + 1 );
			f = fopen( file, "wb" );
			if( f ) {
				if( !GetSwitch( "rawcol", swtc, swtn ) )
				{
					uint8_t* colSrc = color;
					uint8_t* colDst = color;
					for( int b = 0; b < (wc*hc); b += 2 ) {
						*colDst++ = (colSrc[ b ] & 0xf) | (colSrc[ b + 1 ] << 4);
					}
					fwrite( color, wc*hc / 2, 1, f );
				}
				else {
					fwrite( color, wc*hc, 1, f );
				}
				fclose( f );
			}
		}

		free( color );
		free( chars );
		free( screen );
		return 0;
	}

	if (GetSwitch("multisprite", swtc, swtn)) {
		return MultiSprite(args, argn, swtc, swtn);
	}

	if( GetSwitch( "screens", swtc, swtn ) )
	{
		if( argn < 2 ) { printf( "Usage:\nGfx -screens <out> <in-.imap.iscr 1> <in-.imap.iscr 2> <in-.imap.iscr 3>\n" ); return 0; }
		const char* out = args[ 1 ];
		char file[ _MAX_PATH ];
		uint8_t refd[ 256 ];
		const char* imapExt = ".chr";
		const char* iscrExt = ".scr";
		const char* chrExt = ".chr";
		const char* scrExt = ".scr";
		const char* oscrExt = ".scrs";
		int err = 0;
		uint8_t *all = (uint8_t*)calloc( 1, 8 * 256 );
		int numChr = 1;

		for( int i = 2; i < argn; ++i ) {
			uint8_t ecm = 0;
			if( args[i][0] == '*' ) {
				memcpy( file, args[ i ]+1, strlen( args[ i ] )-1 );
				memcpy( file + strlen( args[ i ] )-1, chrExt, sizeof( chrExt ) + 1 );
				ecm = 1;
			} else {
				memcpy( file, args[ i ], strlen( args[ i ] ) );
				memcpy( file + strlen( args[ i ] ), imapExt, sizeof( imapExt ) + 1 );
			}
			FILE* f = fopen( file, "rb" );
			if( f ) {
				fseek( f, 0, SEEK_END );
				size_t size = ftell( f );
				fseek( f, 0, SEEK_SET );
				uint8_t* data = (uint8_t*)malloc( size );
				fread( data, size, 1, f );
				fclose( f );

				size_t added_size = 0;

				if( args[i][0] == '*' ) { memcpy( file + strlen( args[ i ] )-1, scrExt, sizeof( scrExt ) + 1 ); }
				else { memcpy( file + strlen( args[ i ] ), iscrExt, sizeof( iscrExt ) + 1 ); }
				FILE* f = fopen( file, "rb" );
				if( f ) {
					fseek( f, 0, SEEK_END );
					size_t sizeScrn = ftell( f );
					fseek( f, 0, SEEK_SET );
					uint8_t* scrn = (uint8_t*)malloc( sizeScrn );
					uint8_t* scrn2 = (uint8_t*)malloc( sizeScrn );
					fread( scrn, sizeScrn, 1, f );
					fclose( f );

					memset( refd, 0, sizeof( refd ) );

					for( size_t a = 0; a < sizeScrn; ++a ) {
						uint8_t c = scrn[ a ];
						if( ecm ) { c &= 0x3f; }
						if( !refd[ c ] ) {
							uint8_t *cmp = data + (int)(c) * 8;
							uint8_t* chr = all;
							uint8_t found = 0;
							uint8_t index = 0;
							for( int s = 0; s < numChr && !found; ++s ) {
								uint8_t same = 1;
								for( int r = 0; r < 8 && same; ++r ) {
									if( chr[ r ] != cmp[ r ] ) { same = 0; }
								}
								if( same ) {
									found++;
									index = (uint8_t)s;
									break;
								}
								chr += 8;
							}
							if( !found ) {
								index = numChr;
								if( numChr < 256 ) {
									uint8_t *newChr = all + (numChr << 3);
									for( int r = 0; r < 8; ++r ) {
										*newChr++ = *cmp++;
									}
									added_size += 8;
									++numChr;
								}
							}
							for( size_t b = 0; b < sizeScrn; ++b ) {
								if( ecm ) {
									if( (scrn[ b ]&0x3f) == c ) { scrn2[ b ] = index | (scrn[b]&0xc0); }
								} else if( scrn[ b ] == c ) { scrn2[ b ] = index; }
							}
							refd[ c ]++;
						}
					}

					printf( "%s added %d bytes, orig %d bytes\n", args[ i ], (int)added_size, (int)size );

					if( ecm ) {
						memcpy( file, args[ i ]+1, strlen( args[ i ] )-1 );
						memcpy( file + strlen( args[ i ] )-1, oscrExt, sizeof( oscrExt ) + 1 );
					} else {
						memcpy( file, args[ i ], strlen( args[ i ] ) );
						memcpy( file + strlen( args[ i ] ), oscrExt, sizeof( oscrExt ) + 1 );
					}
					FILE* f = fopen( file, "wb" );
					if( f ) {
						fwrite( scrn2, sizeScrn, 1, f );
						fclose( f );
					}
					else {
						printf( "Could not open file %s for writing!\n", file );
						err++;
					}
					free( scrn );
					free( scrn2 );
				}
				free( data );
			}
			else {
				printf( "Could not open file %s for reading!\n", file );
				err++;
			}
		}
		if( !err ) {
			FILE* f = fopen( out, "wb" );
			if( f ) {
				fwrite( all, numChr * 8, 1, f );
				fclose( f );
			}
			else { printf( "couldn't open %s for writing\n", out ); err++; }
			printf( "Total chars in bundle %s = %d\n", out, numChr );
		}
		return err;
	}

	if( GetSwitch( "bundle", swtc, swtn ) )
	{
		if( argn < 2 ) { printf( "Usage:\nGfx -bundle <out> <in set 1> <in set 2> <in set 3>\n"); return 0; }
		const char* out = args[1];
		char file[ _MAX_PATH ];
		const char* extChr = ".chr";
		const char* extChrMap = ".chrmap";
		int err = 0;
		uint8_t *all = (uint8_t*)malloc( 8 * 256 );
		int numChr = 0;

		for( int i = 2; i < argn; ++i ) {
			memcpy( file, args[i], strlen(args[i]) );
			memcpy( file + strlen( args[ i ] ), extChr, sizeof( extChr ) + 1 );
			FILE* f = fopen( file, "rb" );
			if( f ) {
				fseek( f, 0, SEEK_END );
				size_t size = ftell( f );
				fseek( f, 0, SEEK_SET );
				uint8_t* data = (uint8_t*)malloc( size );
				fread( data, size, 1, f );
				fclose( f );

				uint8_t* map = (uint8_t*)malloc( size / 8 );
				uint8_t* cmp = data;
				int matches = 0;
				for( size_t c = 0, n = size / 8; c < n; ++c ) {
					uint8_t* chr = all;
					uint8_t found = 0;
					for( int s=0; s < numChr && !found; ++s ) {
						uint8_t same = 1;
						for( int r = 0; r < 8 && same; ++r ) {
							if( chr[r] != cmp[r] ) { same = 0; }
						}
						if( same ) {
							found++;
							map[ c ] = (uint8_t)s;
							break;
						}
						chr += 8;
					}
					if( !found ) {
						map[ c ] = numChr;
						if( numChr < 256 ) {
							for( int r = 0; r < 8; ++r ) {
								all[ ( numChr<<3 ) + r ] = cmp[ r ];
							}
							++numChr;
						}
					}
					cmp += 8;
				}
				memcpy( file, args[ i ], strlen( args[ i ] ) );
				memcpy( file + strlen( args[ i ] ), extChrMap, sizeof( extChrMap ) + 1 );
				FILE* f = fopen( file, "wb" );
				if( f ) {
					fwrite( map, size/8, 1, f );
					fclose( f );
				} else {
					printf( "Could not open file %s for reading!\n", file );
					err++;
					}
				free( map );
			} else {
				printf( "Could not open file %s for reading!\n", file );
				err++;
			}
		}
		if( !err ) {
			FILE* f = fopen( out, "wb" );
			if( f ) {
				fwrite( all, numChr * 8, 1, f );
				fclose( f );
			} else { printf( "couldn't open %s for writing\n", out ); err++; }
			printf("Total chars in bundle %s = %d\n", out, numChr );
		}
		return err;
	}

	if (GetSwitch("texthires", swtc, swtn)) {
		int w, h;
		uint8_t* img = LoadPicture(args[1], &w, &h);

		int wc = w/8;
		int hc = h/8;

		uint8_t* screen = (uint8_t*)malloc(wc * hc);
		uint8_t* chars = (uint8_t*)malloc(256*8);
		uint8_t* color = (uint8_t*)malloc(wc * hc);
		const char* bgstr = GetSwitch("bg", swtc, swtn);
		uint8_t bg = bgstr ? atoi(bgstr) : 0;

		int nunChars = 1;	// first char is empty!
		memset(chars, 0, 8);

		// prioritize background color 0, then 1 etc.
		for (int y = 0; y < hc; ++y) {
			for (int x = 0; x < wc; ++x) {
				uint8_t* s = img+y*8*w+x*8;
				int fgc = -1;
				uint8_t chr[8] = { 0 };
				for (int yp = 0; yp < 8; ++yp) {
					uint8_t m = 0x80;
					for (int xp = 0; xp < 8; ++xp) {
						int c = s[xp+yp * w];
						if (c!=bg) {
							fgc = c;
							chr[yp] |= m;
						}
						m >>= 1;
					}
				}
				if (fgc>=0) {
					int c = 0;
					for (; c < nunChars; ++c) {
						int same = 1;
						for (int b = 0; b<8; ++b) {
							if (chars[c*8+b]!=chr[b]) {
								same = 0;
								break;
							}
						}
						if (same) { break; }
					}
					if (c==nunChars) {
						if (nunChars < 256) {
							for (int b = 0; b<8; ++b) { chars[c*8+b] = chr[b]; }
							++nunChars;
						} else { c = 0; }
					}
					screen[x+y*wc] = c;
					color[y*wc+x] = fgc;
				} else {
					screen[x+y * wc] = 0;
					color[x+y * wc] = 0;
				}
			}
		}

		printf("Used chars = %d\n", nunChars);

		const char* out = GetSwitch("out", swtc, swtn);
		if (out) {
			size_t outLen = strlen(out);
			char file[_MAX_PATH];
			const char* extChr = ".chr";
			const char* extScr = ".scr";
			const char* extCol = ".col";
			memcpy(file, out, outLen);
			memcpy(file+outLen, extChr, strlen(extChr)+1);
			FILE* f = fopen(file, "wb");
			if (f) {
				if (GetSwitch("skip0", swtc, swtn)) {
					fwrite(chars+8, nunChars*8-8, 1, f);
				} else {
					fwrite(chars, nunChars*8, 1, f);
				}
				fclose(f);
			}

			memcpy(file, out, outLen);
			memcpy(file+outLen, extScr, strlen(extScr)+1);
			f = fopen(file, "wb");
			if (f) {
				fwrite(screen, wc*hc, 1, f);
				fclose(f);
			}

			memcpy(file, out, outLen);
			memcpy(file+outLen, extCol, strlen(extCol)+1);
			f = fopen(file, "wb");
			if (f) {
				if (!GetSwitch("rawcol", swtc, swtn)) {
					uint8_t* colSrc = color;
					uint8_t* colDst = color;
					for (int b = 0; b < (wc*hc); b += 2) {
						*colDst++ = (colSrc[b]&0xf)|(colSrc[b+1]<<4);
					}
					fwrite(color, wc*hc/2, 1, f);
				} else {
					fwrite(color, wc*hc, 1, f);
				}
				fclose(f);
			}
		}

		free(color);
		free(chars);
		free(screen);
		return 0;
	}


	if( argn < 6 ) {
		printf( "Usage:\nGfx img.png bg0 bg1 bg2 bg3\n" );
		return 0;
	}

	uint8_t bg[4] = { (uint8_t)atoi( args[2] ), (uint8_t)atoi( args[3] ), (uint8_t)atoi( args[4] ), (uint8_t)atoi( args[5] ) };

	printf("Converting port %s with background colors: %d, %d, %d, %d (water): ", args[1], bg[0], bg[1], bg[2], bg[3] );

	int w, h;
	uint8_t* img = LoadPicture( args[ 1 ], &w, &h );

	int wc = w / 8;
	int hc = h / 8;

	uint8_t* screen = (uint8_t*)malloc( wc * hc );
	uint8_t* chars = (uint8_t*)malloc( 256 * 8 );
	uint8_t* color = (uint8_t*)malloc( wc * hc );

	int nunChars = 1;	// first char is empty!
	memset( chars, 0, 8 );

	// prioritize background color 0, then 1 etc.
	for( int y = 0; y < hc; ++y ) {
		for( int x = 0; x < wc; ++x ) {
			uint8_t* s = img + y * 8 * w + x * 8;
			int bgc = -1;
			int fgc = -1;
			for( int yp = 0; yp < 8; ++yp ) {
				for( int xp = 0; xp < 8; ++xp ) {
					int c = s[ xp + yp * w ];
					if( c == fgc ) {}
					else if( fgc < 0 ) { fgc = c; }
					else if( bgc < 0 ) { bgc = c; break; }
				}
			}
			int bgi = -1;
			for( int b = 0; b < 4; ++b ) {
				if ( bgc == bg[b] ) { bgi = b; break; }
				else if( fgc == bg[ b ] ) {
					fgc = bgc;
					bgi = b;
					break;
				}
			}
			if (bgi < 0) {
				bgi = 0;
			}
			if( bgi >= 0 ) {
				uint8_t chr[8] = {0};
				for( int yp = 0; yp < 8; ++yp ) {
					for( int xp = 0; xp < 8; ++xp ) {
						if( s[ (xp) + yp * w ] == fgc ) {
							chr[yp] |= 1<<(7-xp);
						}
					}
				}
				int c = 0;
				for( ; c < nunChars; ++c ) {
					int same = 1;
					for( int b=0; b<8; ++b ) {
						if( chars[ c * 8 + b ] != chr[ b ] ) {
							same = 0;
							break;
						}
					}
					if( same ) { break; }
				}
				if( c == nunChars ) {
					if( nunChars < 256 ) {
						for( int b=0; b<8; ++b ) { chars[c*8+b] = chr[b]; }
						++nunChars;
					} else { c = 0; }
				}
				screen[ x + y * wc ] = c + 64 * bgi;
				color[ x + y * wc ] = fgc;
			} else {
				screen[ x + y * wc ] = 0;
				color[ x + y * wc ] = 0;
			}
		}
	}

	printf( "Used chars = %d\n", nunChars );

	if( nunChars < 64 ) {
		//		char* GetSwitch( const char* match, char** swtc, int swtn )
		const char* out = GetSwitch( "out", swtc, swtn );
		if( out ) {
			size_t outLen = strlen( out );
			char file[ _MAX_PATH ];
			const char* extChr = ".chr";
			const char* extScr = ".scr";
			const char* extCol = ".col";
			memcpy( file, out, outLen );
			memcpy( file + outLen, extChr, sizeof( extChr ) + 1 );
			FILE* f = fopen( file, "wb" );
			if( f ) {
				if( GetSwitch( "skip0", swtc, swtn ) ) {
					fwrite( chars + 8, nunChars * 8 - 8, 1, f );
				} else {
					fwrite( chars, nunChars * 8, 1, f );
				}
				fclose( f );
			}

			memcpy( file, out, outLen );
			memcpy( file + outLen, extScr, sizeof( extScr ) + 1 );
			f = fopen( file, "wb" );
			if( f ) {
				fwrite( screen, wc*hc, 1, f );
				fclose( f );
			}

			memcpy( file, out, outLen );
			memcpy( file + outLen, extCol, sizeof( extCol ) + 1 );
			f = fopen( file, "wb" );
			if( f ) {
				if( !GetSwitch( "rawcol", swtc, swtn ) )
				{
					uint8_t* colSrc = color;
					uint8_t* colDst = color;
					for( int b = 0; b < (wc*hc); b += 2 ) {
						*colDst++ = (colSrc[ b ] & 0xf) | (colSrc[ b + 1 ] << 4);
					}
					fwrite( color, wc*hc / 2, 1, f );
				} else {
					fwrite( color, wc*hc, 1, f );
				}
				fclose( f );
			}
		}
	}

	free( color );
	free( chars );
	free( screen );

	return 0;
}

