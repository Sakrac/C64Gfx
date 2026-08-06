#ifndef GPX_LOAD_H
#define GPX_LOAD_H

#include <stddef.h>
#include <stdint.h>

typedef enum gpx_mode {
    GPX_MODE_BITMAP = 0,
    GPX_MODE_MULTICOLOR_BITMAP = 1,
    GPX_MODE_SPRITE = 2,
    GPX_MODE_MULTICOLOR_SPRITE = 3,
    GPX_MODE_CHAR = 4,
    GPX_MODE_MULTICOLOR_CHAR = 5,
    GPX_MODE_UNUSED1 = 6,
    GPX_MODE_UNUSED2 = 7,
    GPX_MODE_UNRESTRICTED = 8,
    GPX_MODE_WIDE_UNRESTRICTED = 9
} gpx_mode;

typedef struct {
    int32_t version;
    int32_t mode;
    int32_t xsize;
    int32_t ysize;
    int32_t mapsize;
    int32_t colorsize;
    int32_t screensize;
    int32_t backbuffers;
    int32_t backbufsel;
    int32_t par;
    int32_t overflow;
	uint8_t background_color;
    uint8_t multicolor0;
    uint8_t multicolor1;
    uint32_t payload_offset;
} gpx_info;

typedef enum {
    GPX_OK = 0,
    GPX_ERR_INVALID_INPUT,
    GPX_ERR_OUT_OF_MEMORY,
    GPX_ERR_DECOMPRESSION_FAILED,
    GPX_ERR_UNSUPPORTED_VERSION,
    GPX_ERR_LEGACY_UNSUPPORTED,
    GPX_ERR_TRUNCATED_METADATA,
    GPX_ERR_BAD_METADATA
} gpx_status;

// Read the the GPX header and metadata, decompressing the payload if necessary. Returns GPX_OK on success, or an error code on failure. The caller is responsible for freeing the out_buffer if the function returns GPX_OK.
gpx_status parse_gpx(const unsigned char *in_data, size_t in_size, gpx_info *info, unsigned char **out_buffer, size_t *out_size);

// Generate a bitmap from the GPX data. Returns GPX_OK on success, or an error code on failure. The caller is responsible for freeing the out_buffer if the function returns GPX_OK.
gpx_status gpx_generate_bitmap(const unsigned char *buffer, size_t buffer_size, const gpx_info *info, unsigned char **out_buffer, size_t *out_size);

// Convert a gpx_status code to a human-readable string.
const char *gpx_status_to_string(gpx_status status);

#endif
