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

This implementation uses the miniz library for decompression
please see miniz/miniz.c for license and other information.

No part of the Pixcen source code is used in this implementation.

This implementation is free to copy, use and modify for any
purpose. It first appeared in https://github.com/sakrac/C64Gfx

*/

#include "gpx_load.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define MINIZ_NO_ARCHIVE_APIS
#define MINIZ_NO_STDIO
#include "miniz/miniz.c"

static int32_t read_le32(const unsigned char *p) {
    return (int32_t)((uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24));
}

static int copy_ascii_string(const unsigned char **p, const unsigned char *end, char *buf, size_t buf_size) {
    size_t i = 0;

    while (*p < end && **p != '\0') {
        if (i + 1 < buf_size) {
            buf[i++] = (char)**p;
        }
        ++(*p);
    }

    if (*p >= end) {
        return 0;
    }

    if (i + 1 < buf_size) {
        buf[i] = '\0';
    }
    ++(*p);
    return 1;
}

static int32_t read_utf16_int(const unsigned char **p, const unsigned char *end) {
    const unsigned char *u = *p;
    int32_t value = 0;
    int sign = 1;

    while (u + 1 < end) {
        uint16_t c = (uint16_t)(u[0] | ((uint16_t)u[1] << 8));
        if (c == 0) {
            break;
        }
        if (c == '-') {
            sign = -1;
        } else if (c >= '0' && c <= '9') {
            value = value * 10 + (c - '0');
        } else {
            break;
        }
        u += 2;
    }

    if (u + 1 <= end) {
        uint16_t c = (uint16_t)(u[0] | ((uint16_t)u[1] << 8));
        if (c == 0) {
            u += 2;
        }
    }

    *p = u;
    return sign * value;
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
            char key[64];
            if (p >= end) {
                free(decompressed);
                return GPX_ERR_TRUNCATED_METADATA;
            }
            if (!copy_ascii_string(&p, end, key, sizeof(key))) {
                free(decompressed);
                return GPX_ERR_TRUNCATED_METADATA;
            }
            if (p + 1 > end) {
                free(decompressed);
                return GPX_ERR_TRUNCATED_METADATA;
            }
            if (strcmp(key, "xsize") == 0) {
                info->xsize = read_utf16_int(&p, end);
            } else if (strcmp(key, "ysize") == 0) {
                info->ysize = read_utf16_int(&p, end);
            } else if (strcmp(key, "mapsize") == 0) {
                info->mapsize = read_utf16_int(&p, end);
            } else if (strcmp(key, "colorsize") == 0) {
                info->colorsize = read_utf16_int(&p, end);
            } else if (strcmp(key, "screensize") == 0) {
                info->screensize = read_utf16_int(&p, end);
            } else if (strcmp(key, "backbuffers") == 0) {
                info->backbuffers = read_utf16_int(&p, end);
            } else if (strcmp(key, "backbufsel") == 0) {
                info->backbufsel = read_utf16_int(&p, end);
            } else if (strcmp(key, "par") == 0) {
                info->par = read_utf16_int(&p, end);
            } else if (strcmp(key, "overflow") == 0) {
                info->overflow = read_utf16_int(&p, end);
            } else {
                read_utf16_int(&p, end);
            }
        }
    }

    if (info->mapsize <= 0 || info->backbuffers <= 0) {
        free(decompressed);
        return GPX_ERR_BAD_METADATA;
    }

    info->payload_offset = (size_t)(p - decompressed);

	if (info->payload_offset + 64 <= (size_t)decompressed_size) {
		size_t slot_size = 64 + (size_t)info->mapsize + (size_t)info->screensize + (size_t)info->colorsize;
		size_t slot_index = (size_t)(info->backbufsel >= 0 ? info->backbufsel : 0);
		size_t slot_offset = 64 + slot_size * slot_index;
		size_t slot_base = info->payload_offset + slot_offset;

		if (slot_base + 64 <= (size_t)decompressed_size) {
			info->background_color = decompressed[slot_base + 60];
			info->char_multicolor0 = decompressed[slot_base + 61];
			info->char_multicolor1 = decompressed[slot_base + 62];
			info->sprite_multicolor0 = decompressed[slot_base + 63];
			info->sprite_multicolor1 = 0;
		}
	}

    *out_buffer = decompressed;
    *out_size = decompressed_size;
    return GPX_OK;
}

uint8_t* create_image(gpx_mode mode, int w, int h, uint8_t *screen, uint8_t *colors, uint8_t *chars)
{
    if (mode == GPX_MODE_SPRITE || mode == GPX_MODE_MULTICOLOR_SPRITE) {
        // For sprite modes, we need to handle differently, but for now, let's just return NULL
        return NULL;
    }

	bool mc = (mode == GPX_MODE_MULTICOLOR_BITMAP || mode == GPX_MODE_MULTICOLOR_CHAR);

    uint8_t* image = (uint8_t*)malloc(w * h * (mc ? 2 : 1));
    if (!image) {
        return NULL;
    }

    if (mc) { w <<= 1; }

    // character based modes
    int cw = w>>3, ch = h>>3;

    for(int y=0; y<h; y++)
    {
        for(int x=0; x<w; x++)
        {
            int index = y * w + x;
            uint8_t pixel_value = 0;
            int cx = x>>3, cy = y>>3;
            int map_offs = cx + cy * cw;

            uint8_t screen_value = screen[map_offs];
            uint8_t color_value = colors[map_offs];
            uint8_t color_map[4] = {0, color_value, 1, 1 };
            uint8_t *chardata = NULL;

            switch(mode)
            {
                case GPX_MODE_MULTICOLOR_BITMAP:
					chardata = chars + 8 * map_offs;
                    color_map[0] = 9;// screen_value & 0x0f;
					color_map[1] = screen_value >> 4;
					color_map[2] = screen_value & 0xf;
					color_map[3] = color_value & 0xf;
					break;
                case GPX_MODE_BITMAP:
                    chardata = chars + 8 * map_offs;
                    color_map[0] = screen_value & 0x0f;
                    color_map[1] = (screen_value >> 4) & 0x0f;
                    break;
                case GPX_MODE_CHAR:
                case GPX_MODE_MULTICOLOR_CHAR:
                    chardata = chars + 8 * screen_value;
                    color_map[1] = color_value & 0x0f;
                    break;
            }
            if (chardata) {
                uint8_t b = chardata[y&7];
                if (mc) {
                    pixel_value = color_map[(b >> ((~x)&6)) & 3];
                } else {
                    pixel_value = color_map[(b >> ((~x)&7)) & 1];
                }
                image[index] = pixel_value;
            }
        }
    }
    return image;
}

gpx_status gpx_build_indexed_bitmap(const unsigned char *buffer, size_t buffer_size, const gpx_info *info, unsigned char **out_buffer, size_t *out_size) {
    size_t decoded_count;
    size_t slot_size;
    size_t slot_index;
    size_t slot_offset;
    size_t map_offset;
    size_t i;
    unsigned char *bitmap;

    if (!buffer || !info || !out_buffer || !out_size) {
        return GPX_ERR_INVALID_INPUT;
    }

    if (info->xsize <= 0 || info->ysize <= 0 || info->mapsize <= 0) {
        return GPX_ERR_BAD_METADATA;
    }

    if (info->backbuffers <= 0) {
        return GPX_ERR_BAD_METADATA;
    }

    slot_size = 64 + (size_t)info->mapsize + (size_t)info->screensize + (size_t)info->colorsize;
    slot_index = (size_t)(info->backbufsel >= 0 ? info->backbufsel : 0);
    if (slot_index >= (size_t)info->backbuffers) {
        return GPX_ERR_BAD_METADATA;
    }
    /* GPX payload layout:
    - a 64-byte prefix before the first slot
    - then one slot per backbuffer, each of size 64 + mapsize + screensize + colorsize
    - inside each slot, the map data starts at offset 64 from the slot base
    */
    slot_offset = 64 + slot_size * slot_index;
    map_offset = info->payload_offset + slot_offset + 64;

    if (map_offset > buffer_size || (size_t)info->mapsize > buffer_size - map_offset) {
        return GPX_ERR_BAD_METADATA;
    }

    const uint8_t* pChars = buffer + map_offset;
    const uint8_t* pColors = pChars + info->mapsize;
    const uint8_t* pScreen = pColors + info->colorsize;

/*
then the map bytes
then the color bytes
then the screen bytes
*/
#if 1
//    uint8_t *screen_data = buffer + info->payload_offset + 64;
//    uint8_t *colors_data = screen_data + info->screensize;
//    uint8_t *chars_data = colors_data + info->colorsize;

    bitmap = create_image((gpx_mode)info->mode, info->xsize, info->ysize, pScreen, pColors, pChars);
	*out_buffer = bitmap;
	*out_size = info->xsize * info->ysize;
	return GPX_OK;
#else

    switch (info->mode) {
        case 0: /* bitmap */
        case 1: /* multicolor bitmap */
        case 8: /* unrestricted */
        case 9: /* wide unrestricted */
            break;
        default:
            return GPX_ERR_BAD_METADATA;
    }

    decoded_count = (size_t)info->xsize * (size_t)info->ysize;
    bitmap = (unsigned char *)malloc(decoded_count);
    if (!bitmap) {
        return GPX_ERR_OUT_OF_MEMORY;
    }

    bool mc = (info->mode == GPX_MODE_MULTICOLOR_BITMAP || info->mode == GPX_MODE_MULTICOLOR_CHAR || info->mode == GPX_MODE_MULTICOLOR_SPRITE);



    for (i = 0; i < decoded_count; ++i) {
        size_t x = i % (size_t)info->xsize;
        size_t y = i / (size_t)info->xsize;
        size_t cell_x;
        size_t cell_y;
        size_t byte_index;
        unsigned char map_byte;

        if (mc) {


            cell_x = x / 4;
            cell_y = y / 8;
			int screen_index = cell_y * (info->xsize>>2) + cell_x;
            byte_index = cell_y * (info->xsize>>2) + cell_x * 8 + (y&7);
            map_byte = buffer[map_offset + byte_index];
			uint8_t screen_byte = pScreen[screen_index];
			uint8_t color_byte = pColors[screen_index];
            uint8_t color_index = (map_byte >> (6 - 2*(x&3))) & 3;
            switch (color_index) {
                case 0:
                    bitmap[i] = 0; // Background color
                    break;
                case 1:
                    bitmap[i] = screen_byte >> 4;
                    break;
                case 2:
                    bitmap[i] = screen_byte & 0xf;
                    break;
                case 3:
                    bitmap[i] = color_byte & 0xf;
                    break;
			}
            bitmap[i] = (unsigned char)((map_byte >> (7 - (x % 8))) & 0x01);
        } else if (info->mode != GPX_MODE_UNRESTRICTED && info->mode != GPX_MODE_WIDE_UNRESTRICTED) {
            cell_x = x / 4;
            cell_y = y / 8;
            byte_index = cell_y * (size_t)info->xsize * 2 + cell_x * 8 + (y % 8);
            map_byte = buffer[map_offset + byte_index];
            bitmap[i] = (unsigned char)((map_byte >> (2 * (3 - (x % 4)))) & 0x03);
        } else {
            bitmap[i] = (unsigned char)(buffer[map_offset + y * (size_t)info->xsize + x] & 0x0F);
        }
    }
	*out_buffer = bitmap;
	*out_size = decoded_count;
	return GPX_OK;
#endif
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
