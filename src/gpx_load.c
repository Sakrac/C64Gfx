/*

This is a re-implementation of the GPX file format parser,
which is used to read and extract metadata from GPX files.
The parser decompresses the input data using the miniz library
and reads various metadata fields such as version, mode,
xsize, ysize, mapsize, colorsize, screensize, backbuffers,
backbufsel, par, and overflow.
The extracted information is stored in a gpx_info structure.

GPX is the proprietary graphics format used by Pixcen that
is developed by John Hammarberg: https://github.com/Hammarberg/pixcen/

Pixcen downloads: https://hammarberg.github.io/pixcen/

This implementation uses the miniz library for decompression,
see miniz/miniz.c for license and other information.

No part of the Pixcen source code is used in this implementation.

This implementation is free to copy, use and modify for any
purpose. It first appeared in https://github.com/sakrac/C64Gfx
(by Carl-Henrik Skårman and is released under the MIT license.)

*/

/* Keep MSVC from warning on fopen in the file writer. */
#define _CRT_SECURE_NO_WARNINGS

#include "gpx_load.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MINIZ_NO_ARCHIVE_APIS
#define MINIZ_NO_STDIO
#include "miniz/miniz.c"

// header + 128 bytes of metadata
#define GPX_MAX_HEADER_SIZE 1024

static const char *gpx_metadata_keys[] = {
	"xsize",
	"ysize",
	"mapsize",
	"colorsize",
	"screensize",
	"backbuffers",
	"backbufsel",
	"par",
	"overflow"
};

uint32_t gpx_assign_u32(uint32_t value, uint8_t *out_bytes) {
	out_bytes[0] = (uint8_t)(value & 0xffu);
	out_bytes[1] = (uint8_t)((value >> 8) & 0xffu);
	out_bytes[2] = (uint8_t)((value >> 16) & 0xffu);
	out_bytes[3] = (uint8_t)((value >> 24) & 0xffu);
	return 4;
}

uint32_t gpx_assign_wide_num(uint32_t value, uint8_t *out_bytes) {
	uint32_t digit_count = 1, digit_value = 1;
	while (value >= (10 * digit_value)) {
		digit_count++;
		digit_value *= 10;
	}
	for (uint32_t index = 0; index < digit_count; ++index) {
		uint32_t digit = (value / digit_value) % 10;
		out_bytes[index * 2] = (uint8_t)('0' + digit);
		out_bytes[index * 2 + 1] = 0;
		digit_value /= 10;
	}
	out_bytes[digit_count * 2] = 0x00;
	out_bytes[digit_count * 2 + 1] = 0x00;
	return ( digit_count + 1 ) * 2;
}

uint32_t gpx_assign_string_wide_num(const char *key, uint32_t value, uint8_t *out_bytes) {
	uint32_t key_length = (uint32_t)strlen(key);
	memcpy(out_bytes, key, key_length);
	out_bytes[key_length] = 0;
	return gpx_assign_wide_num(value, out_bytes + key_length + 1) + key_length + 1;
}

uint32_t gpx_make_header(uint8_t *header, const gpx_info *info) {
	uint32_t offset = 0;
	offset += gpx_assign_u32((uint32_t)info->version, header + offset);
	offset += gpx_assign_u32((uint32_t)info->mode, header + offset);
	offset += gpx_assign_u32(9, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[0], (uint32_t)info->xsize, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[1], (uint32_t)info->ysize, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[2], (uint32_t)info->mapsize, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[3], (uint32_t)info->colorsize, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[4], (uint32_t)info->screensize, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[5], (uint32_t)info->backbuffers, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[6], (uint32_t)info->backbufsel, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[7], (uint32_t)info->par, header + offset);
	offset += gpx_assign_string_wide_num(gpx_metadata_keys[8], (uint32_t)info->overflow, header + offset);
	return offset;
}

uint8_t* gpx_create(const GpxData *data, const gpx_info *info, size_t *out_size) {

	size_t total_size = GPX_MAX_HEADER_SIZE + info->mapsize + info->screensize + info->colorsize;

	uint8_t *out = (uint8_t *)calloc(1, total_size);

	uint32_t offset = gpx_make_header(out, info);

	offset += 64; // 0 header
	out[offset + 60] = info->background_color;
	out[offset + 61] = info->multicolor0;
	out[offset + 62] = info->multicolor1;

	offset += 64; // 1 header

	if (info->mapsize) { memcpy(out + offset, data->pChars, info->mapsize); }
	offset += info->mapsize;
	if (info->colorsize) { memcpy(out + offset, data->pColors, info->colorsize); }
	offset += info->colorsize;
	if (info->screensize) { memcpy(out + offset, data->pScreen, info->screensize); }
	offset += info->screensize;

	mz_ulong compressed_size = compressBound((mz_ulong)offset);
	uint8_t *compressed_data = (uint8_t *)malloc(compressed_size);
	if (!compressed_data) {
		free(out);
		return NULL;
	}

	compressed_size = (mz_ulong)compressed_size;
	uint32_t status = compress2(compressed_data, &compressed_size, out, (mz_ulong)offset, MZ_BEST_COMPRESSION);
	if (status != MZ_OK) {
		free(out);
		free(compressed_data);
		return NULL;
	}

	*out_size = compressed_size;
	free(out);
	return compressed_data;
}

static inline int32_t read_le32(const unsigned char *p) {
    return (int32_t)((uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24));
}

gpx_status parse_gpx(const unsigned char *in_data, size_t in_size, gpx_info *info, unsigned char **out_buffer, size_t *out_size) {
	const unsigned char *p;
	const unsigned char *end;
	unsigned char *decompressed = NULL;
	mz_ulong decompressed_size = 0;
	int status;

	if (!in_data || !info || !out_buffer || !out_size || in_size < 16) {
		return GPX_ERR_INVALID_INPUT;
	}

	decompressed_size = (mz_ulong)(in_size * 8);
	decompressed = (unsigned char *)malloc(decompressed_size);
	if (!decompressed) {
		return GPX_ERR_OUT_OF_MEMORY;
	}

	for (;;) {
		mz_ulong tmp_len = decompressed_size;
		status = mz_uncompress(decompressed, &tmp_len, in_data, (mz_ulong)in_size);
		if (status == MZ_OK) {
			decompressed_size = tmp_len;
			break;
		}
		if (status == MZ_BUF_ERROR) {
			decompressed_size *= 2;
			decompressed = (unsigned char *)realloc(decompressed, decompressed_size);
			if (!decompressed) {
				free(decompressed);
				return GPX_ERR_OUT_OF_MEMORY;
			}
			continue;
		}
		free(decompressed);
		return GPX_ERR_DECOMPRESSION_FAILED;
	}

	p = decompressed;
	end = decompressed + decompressed_size;

	if (decompressed_size < 8) {
		free(decompressed);
		return GPX_ERR_INVALID_INPUT;
	}

	memset(info, 0, sizeof(*info));
	info->version = read_le32(p);
	p += 4;
	info->mode = read_le32(p);
	p += 4;

	if (info->version > 4) {
		free(decompressed);
		return GPX_ERR_UNSUPPORTED_VERSION;
	}

	if (info->version == 0) {
		free(decompressed);
		return GPX_ERR_LEGACY_UNSUPPORTED;
	}

	if (info->version >= 4) {
		int32_t metacount = read_le32(p);
		int32_t i;
		p += 4;

		for (i = 0; i < metacount; ++i) {
			if (p >= end) {
				free(decompressed);
				return GPX_ERR_TRUNCATED_METADATA;
			}
			const char* key = (const char*)p;
			while (p < end && *p) { p++; }
			if(p==end) {
				free(decompressed);
				return GPX_ERR_TRUNCATED_METADATA;
			}
			p++; // skip 0

			int32_t arg = 0;
			while(p<end && *p) {
				if (*p < '0' || *p>'9') {
					free(decompressed);
					return GPX_ERR_TRUNCATED_METADATA;
				}
				arg = arg * 10 + (*p - '0');
				p += 2;
			}
			if (p >= (end-1)) {
				free(decompressed);
				return GPX_ERR_TRUNCATED_METADATA;
			}
			p += 2;

			if (strcmp(key, "xsize") == 0) {
				info->xsize = arg;
			} else if (strcmp(key, "ysize") == 0) {
				info->ysize = arg;
			} else if (strcmp(key, "mapsize") == 0) {
				info->mapsize = arg;
			} else if (strcmp(key, "colorsize") == 0) {
				info->colorsize = arg;
			} else if (strcmp(key, "screensize") == 0) {
				info->screensize = arg;
			} else if (strcmp(key, "backbuffers") == 0) {
				info->backbuffers = arg;
			} else if (strcmp(key, "backbufsel") == 0) {
				info->backbufsel = arg;
			} else if (strcmp(key, "par") == 0) {
				info->par = arg;
			} else if (strcmp(key, "overflow") == 0) {
				info->overflow = arg;
			} else if (strcmp(key, "autoselect") == 0) {
				info->autoselect = (uint8_t)arg;
			} else if (strcmp(key, "cellsnap") == 0) {
				info->cellsnap = (uint8_t)arg;
			} else if (strcmp(key, "palette") == 0) {
				info->palette = arg;
			} else if (strcmp(key, "pcol") == 0) {
				info->pcol = arg;
			} else if (strcmp(key, "pmap") == 0) {
				info->pmap = arg;
			} else if (strcmp(key, "pscr") == 0) {
				info->pscr = arg;
#ifdef _DEBUG
			} else {
				printf("GPX Key: %s, Arg: %d\n", key, arg);
#endif
			}
		}
	}

	// give multicolor double width so it matches Pixcen exported png
	switch(info->mode) {
		case GPX_MODE_MULTICOLOR_BITMAP:
		case GPX_MODE_MULTICOLOR_SPRITE:
		case GPX_MODE_MULTICOLOR_CHAR:
		case GPX_MODE_WIDE_UNRESTRICTED:
			info->xsize *= 2;
			break;
		default:
			break;
	}

	if (info->mapsize <= 0 || info->backbuffers <= 0) {
		free(decompressed);
		return GPX_ERR_BAD_METADATA;
	}

	info->payload_offset = (uint32_t)(p - decompressed);

	if (info->payload_offset + 64 <= (size_t)decompressed_size) {
		size_t slot_size = 64 + (size_t)info->mapsize + (size_t)info->screensize + (size_t)info->colorsize;
		size_t slot_index = (size_t)(info->backbufsel >= 0 ? info->backbufsel : 0);
		size_t slot_offset = 64 + slot_size * slot_index;
		size_t slot_base = info->payload_offset + slot_offset;

		if (slot_base + 64 <= (size_t)decompressed_size) {
			info->background_color = decompressed[slot_base + 60];
			info->multicolor0 = decompressed[slot_base + 61];
			info->multicolor1 = decompressed[slot_base + 62];
		}
	}

	*out_buffer = decompressed;
	*out_size = decompressed_size;
	return GPX_OK;
}

uint8_t* create_image(const GpxData data, size_t *out_size) {
	const uint8_t* chars = data.pChars;
	const uint8_t* colors = data.pColors;
	const uint8_t* screen = data.pScreen;
	int w = data.width;
	int h = data.height;
	gpx_mode mode = data.mode;
	uint8_t bg = data.background_color;
	uint8_t mc0 = data.multicolor0;
	uint8_t mc1 = data.multicolor1;

	uint8_t* image = (uint8_t*)malloc(w * h);
	if (!image) {
		return NULL;
	}
	if (out_size) { *out_size = (size_t)(w * h); }

	// Unrestricted modes
	if (mode == GPX_MODE_UNRESTRICTED)
	{
		memcpy(image, chars, w * h);
		return image;
	} else if (mode == GPX_MODE_WIDE_UNRESTRICTED) {
		for(int i=0, n=w*h; i<n; i++) {
			image[i] = chars[i>>1];
		}
		return image;
	}

	bool mc = mode == GPX_MODE_MULTICOLOR_BITMAP || mode == GPX_MODE_MULTICOLOR_CHAR || mode == GPX_MODE_MULTICOLOR_SPRITE;

	// Sprite based modes
	if (mode == GPX_MODE_SPRITE || mode == GPX_MODE_MULTICOLOR_SPRITE) {
		uint8_t color_lookup[4] = { bg, mc0, 0, mc1 };
		int sw = w / 24;
		for (int sy = 0, sh = h / 21; sy < sh; sy++) {
			for (int sx = 0; sx < sw; sx++) {
				color_lookup[mc ? 2 : 1] = *colors++;
				for (int y = 0; y < 21; y++) {
					for (int xb = 0; xb < 3; xb++) {
						uint8_t b = *chars++;
						for (int x = 0; x < 8; x++) {
							uint8_t color_index = mc ?
								((b >> (~x & 6)) & 3) :
								((b >> (~x & 7)) & 1);
							image[(sy * 21 + y) * w + (sx * 24 + xb * 8 + x)] = color_lookup[color_index];
						}
					}
				}
				chars++; // 1 extra byte per sprite
			}
		}
		return image;
	}

	// character based modes
	int cw = w>>3;
	for (int cy = 0, ch = h>>3; cy < ch; cy++) {
		for (int cx = 0; cx < cw; cx++) {
			const uint8_t screen_value = *screen++;
			const uint8_t color_value = *colors++;
			const uint8_t* chardata = chars;
			chars += 8;

			uint8_t color_map[4] = { bg, color_value, mc0, mc1 };
			switch (mode) {
				case GPX_MODE_MULTICOLOR_BITMAP:
					color_map[1] = screen_value >> 4;
					color_map[2] = screen_value & 0xf;
					color_map[3] = color_value & 0xf;
					break;
				case GPX_MODE_BITMAP:
					color_map[0] = screen_value & 0x0f;
					color_map[1] = (screen_value >> 4) & 0x0f;
					break;
				case GPX_MODE_CHAR:
					color_map[1] = color_value & 0x0f;
					break;
				case GPX_MODE_MULTICOLOR_CHAR:
					color_map[1] = mc0;
					color_map[2] = mc1;
					color_map[3] = color_value & 0xf;
					break;
				default:
					break;
			}
			for (int y = 0; y < 8; y++) {
				int index = (cy * 8 + y) * w + (cx * 8);
				for (int x = 0; x < 8; x++) {
					uint8_t b = chardata[y];
					if (mc) {
						image[index] = color_map[(b >> ((~x) & 6)) & 3];
					} else {
						image[index] = color_map[(b >> ((~x) & 7)) & 1];
					}
					++index;
				}
			}
		}
	}
	return image;
}

GpxData gpx_extract_data(const unsigned char *buffer, size_t buffer_size, const gpx_info *info) {
	GpxData data = {0};
	size_t slot_size;
	size_t slot_index;
	size_t slot_offset;
	size_t map_offset;
	if (!buffer || !info) {
		return data;
	}
	if (info->xsize <= 0 || info->ysize <= 0 || info->mapsize <= 0) {
		return data;
	}
	if (info->backbuffers <= 0) {
		return data;
	}
	slot_size = 64 + (size_t)info->mapsize + (size_t)info->screensize + (size_t)info->colorsize;
	slot_index = (size_t)(info->backbufsel >= 0 ? info->backbufsel : 0);
	if (slot_index >= (size_t)info->backbuffers) {
		return data;
	}
	slot_offset = 64 + slot_size * slot_index;
	map_offset = info->payload_offset + slot_offset + 64;
	if (map_offset > buffer_size || (size_t)info->mapsize > buffer_size - map_offset) {
		return data;
	}
	data.pChars = buffer + map_offset;
	data.pColors = data.pChars + info->mapsize;
	data.pScreen = data.pColors + info->colorsize;
	data.width = (uint16_t)info->xsize;
	data.height = (uint16_t)info->ysize;
	data.mode = (gpx_mode)info->mode;
	data.background_color = info->background_color;
	data.multicolor0 = info->multicolor0;
	data.multicolor1 = info->multicolor1;
	return data;
}

gpx_status gpx_generate_bitmap(const unsigned char *buffer, size_t buffer_size, const gpx_info *info, unsigned char **out_buffer, size_t *out_size) {

	*out_buffer = create_image(gpx_extract_data(buffer, buffer_size, info), out_size);
	return out_buffer ? GPX_OK : GPX_ERR_OUT_OF_MEMORY;
}

const char *gpx_status_to_string(gpx_status status) {
	switch (status) {
		case GPX_OK: return "ok";
		case GPX_ERR_INVALID_INPUT: return "invalid input";
		case GPX_ERR_OUT_OF_MEMORY: return "out of memory";
		case GPX_ERR_DECOMPRESSION_FAILED: return "decompression failed";
		case GPX_ERR_UNSUPPORTED_VERSION: return "unsupported version";
		case GPX_ERR_LEGACY_UNSUPPORTED: return "legacy GPX unsupported";
		case GPX_ERR_TRUNCATED_METADATA: return "truncated metadata";
		case GPX_ERR_BAD_METADATA: return "bad metadata";
		default: return "unknown error";
	}
}
