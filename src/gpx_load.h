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

typedef struct gpx_data {
	const uint8_t* pChars;
	const uint8_t* pColors;
	const uint8_t* pScreen;
	int16_t width;
	int16_t height;
	gpx_mode mode;
	uint8_t background_color;
	uint8_t multicolor0;
	uint8_t multicolor1;
} gpx_data;

// Read the the GPX header and metadata, Returns GPX_OK on success, or an error code on failure. The caller is responsible for freeing *data if the function returns GPX_OK.
gpx_status parse_gpx(const unsigned char* in_data, size_t in_size, gpx_data** data);

// Generate a bitmap from the GPX data. Returns GPX_OK on success, or an error code on failure. The caller is responsible for freeing the out_buffer if the function returns GPX_OK.
gpx_status gpx_generate_bitmap(const gpx_data* data, unsigned char** out_buffer, size_t* out_size);

// Create a Pixcen-compatible GPX file in memory
uint8_t* gpx_create(const gpx_data* data, size_t* out_size);

// Convert a gpx_status code to a human-readable string.
const char* gpx_status_to_string(gpx_status status);

#endif
