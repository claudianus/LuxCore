// cyCodeBase by Cem Yuksel
// [www.cemyuksel.com]
//-------------------------------------------------------------------------------
///
/// \file		cyHairFile.h 
/// \author		Cem Yuksel
/// \version	1.0
/// \date		April 4, 2007
///
/// \brief class for HAIR file type
///
/// Extended with a couple of additional features for LuxRender project.
///
//-------------------------------------------------------------------------------

// NOTE: this file is included in LuxCore so any external dependency must be
// avoided here

#ifndef _CY_HAIR_FILE_H_INCLUDED_
#define _CY_HAIR_FILE_H_INCLUDED_

//-------------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <luxrays/utils/exportdefs.h>

namespace luxrays {

//-------------------------------------------------------------------------------

#define CY_HAIR_FILE_SEGMENTS_BIT		1
#define CY_HAIR_FILE_POINTS_BIT			2
#define CY_HAIR_FILE_THICKNESS_BIT		4
#define CY_HAIR_FILE_TRANSPARENCY_BIT	8
#define CY_HAIR_FILE_COLORS_BIT			16
#define CY_HAIR_FILE_UVS_BIT			32

#define CY_HAIR_FILE_INFO_SIZE			88

// File read errors
#define CY_HAIR_FILE_ERROR_CANT_OPEN_FILE			-1
#define CY_HAIR_FILE_ERROR_CANT_READ_HEADER			-2
#define	CY_HAIR_FILE_ERROR_WRONG_SIGNATURE			-3
#define	CY_HAIR_FILE_ERROR_READING_SEGMENTS			-4
#define	CY_HAIR_FILE_ERROR_READING_POINTS			-5
#define	CY_HAIR_FILE_ERROR_READING_THICKNESS		-6
#define	CY_HAIR_FILE_ERROR_READING_TRANSPARENCY		-7
#define	CY_HAIR_FILE_ERROR_READING_COLORS			-8
#define	CY_HAIR_FILE_ERROR_READING_UVS				-9

//-------------------------------------------------------------------------------

/// Hair file header
struct cyHairFileHeader
{
	char			signature[4];	///< This should be "HAIR"
	unsigned int	hair_count;		///< number of hair strands
	unsigned int	point_count;	///< total number of points of all strands
	unsigned int	arrays;			///< bit array of data in the file

	unsigned int	d_segments;		///< default number of segments of each strand
	float			d_thickness;	///< default thickness of hair strands
	float			d_transparency;	///< default transparency of hair strands
	float			d_color[3];		///< default color of hair strands

	std::array<char, CY_HAIR_FILE_INFO_SIZE>	info;	///< information about the file
};

//-------------------------------------------------------------------------------

/// HAIR file class
CPP_EXPORT class CPP_API cyHairFile
{
public:
	cyHairFile() { Initialize(); }
	~cyHairFile() { Initialize(); }


	//////////////////////////////////////////////////////////////////////////
	///@name Constant Data Access Methods
	
	const cyHairFileHeader& GetHeader() const { return header; }				///< Use this method to access header data.
	std::span<const unsigned short> GetSegmentsArray() const { return segments; }	///< Returns segments array (segment count for each hair strand).
	std::span<const float> GetPointsArray() const { return points; }				///< Returns points array (xyz coordinates of each hair point).
	std::span<const float> GetThicknessArray() const { return thickness; }			///< Returns thickness array (thickness at each hair point}.
	std::span<const float> GetTransparencyArray() const { return transparency; }	///< Returns transparency array (transparency at each hair point).
	std::span<const float> GetColorsArray() const { return colors; }				///< Returns colors array (rgb color at each hair point).
	std::span<const float> GetUVsArray() const { return uvs; }						///< Returns uvs array (uv at each hair point).


	//////////////////////////////////////////////////////////////////////////
	///@name Data Access Methods

	unsigned short* GetSegmentsArray() { return segments.data(); }				///< Returns segments array (segment count for each hair strand).
	float* GetPointsArray() { return points.data(); }							///< Returns points array (xyz coordinates of each hair point).
	float* GetThicknessArray() { return thickness.data(); }						///< Returns thickness array (thickness at each hair point}.
	float* GetTransparencyArray() { return transparency.data(); }				///< Returns transparency array (transparency at each hair point).
	float* GetColorsArray() { return colors.data(); }							///< Returns colors array (rgb color at each hair point).
	float* GetUVsArray() { return uvs.data(); }									///< Returns uvs array (uv at each hair point).


	//////////////////////////////////////////////////////////////////////////
	///@name Methods for Setting Array Sizes
	
	/// Deletes all arrays and initializes the header data.
	void Initialize()
	{
		segments.clear();
		points.clear();
		colors.clear();
		thickness.clear();
		transparency.clear();
		uvs.clear();

		header.signature[0] = 'H';
		header.signature[1] = 'A';
		header.signature[2] = 'I';
		header.signature[3] = 'R';
		header.hair_count = 0;
		header.point_count = 0;
		header.arrays = 0;	// no arrays
		header.d_segments = 0;
		header.d_thickness = 1.0f;
		header.d_transparency = 0.0f;
		header.d_color[0] = 1.0f;
		header.d_color[1] = 1.0f;
		header.d_color[2] = 1.0f;
		std::ranges::fill(std::span(header.info), '\0');
	}

	/// Sets the hair count, re-allocates segments array if necessary.
	void SetHairCount( int count )
	{
		header.hair_count = count;
		if ( not segments.empty() ) {
			segments.clear();
			segments.resize(header.hair_count);
		}
	}

	// Sets the point count, re-allocates points, thickness, transparency, and colors arrays if necessary.
	void SetPointCount(size_t count)
	{
		header.point_count = count;

		auto reset = [&](std::vector<float>& arr, size_t dimension) {
			if (arr.empty()) return;
			arr.clear();
			arr.resize(dimension * count);
		};

		reset(points, 3);
		reset(thickness, 1);
		reset(transparency, 1);
		reset(colors, 3);
		reset(uvs, 2);
	}

	/// Use this function to allocate/delete arrays.
	/// Before you call this method set hair count and point count.
	/// Note that a valid HAIR file should always have points array.
	void SetArrays( int array_types )
	{
		header.arrays = array_types;
		if ( header.arrays & CY_HAIR_FILE_SEGMENTS_BIT && segments.empty() ) segments.resize(header.hair_count);
		if ( ! (header.arrays & CY_HAIR_FILE_SEGMENTS_BIT) && !segments.empty() ) segments.clear();

		if ( header.arrays & CY_HAIR_FILE_POINTS_BIT && points.empty() ) points.resize(header.point_count * 3);
		if ( ! (header.arrays & CY_HAIR_FILE_POINTS_BIT) && !points.empty() ) points.clear();

		if ( header.arrays & CY_HAIR_FILE_THICKNESS_BIT && thickness.empty() ) thickness.resize(header.point_count);
		if ( ! (header.arrays & CY_HAIR_FILE_THICKNESS_BIT) && !thickness.empty() ) thickness.clear();

		if ( header.arrays & CY_HAIR_FILE_TRANSPARENCY_BIT && transparency.empty() ) transparency.resize(header.point_count);
		if ( ! (header.arrays & CY_HAIR_FILE_TRANSPARENCY_BIT) && !transparency.empty() ) transparency.clear();

		if ( header.arrays & CY_HAIR_FILE_COLORS_BIT && colors.empty() ) colors.resize(header.point_count*3);
		if ( ! (header.arrays & CY_HAIR_FILE_COLORS_BIT) && !colors.empty() ) colors.clear();

		if ( header.arrays & CY_HAIR_FILE_UVS_BIT && uvs.empty() ) uvs.resize(header.point_count*2);
		if ( ! (header.arrays & CY_HAIR_FILE_UVS_BIT) && !uvs.empty() ) uvs.clear();

	}

	/// Sets default number of segments for all hair strands, which is used if segments array does not exist.
	void SetDefaultSegmentCount( int s ) { header.d_segments = s; }

	/// Sets default hair strand thickness, which is used if thickness array does not exist.
	void SetDefaultThickness( float t ) { header.d_thickness = t; }

	/// Sets default hair strand transparency, which is used if transparency array does not exist.
	void SetDefaultTransparency( float t ) { header.d_transparency = t; }

	/// Sets default hair color, which is used if color array does not exist.
	void SetDefaultColor( float r, float g, float b ) { header.d_color[0]=r; header.d_color[1]=g; header.d_color[2]=b; }

	//////////////////////////////////////////////////////////////////////////
	///@name Load and Save Methods

	/// Loads hair data from the given HAIR file.
	int LoadFromFile( const char *filename )
	{
		Initialize();

		FILE *fp;
		fp = fopen( filename, "rb" );
		if ( fp == NULL ) return CY_HAIR_FILE_ERROR_CANT_OPEN_FILE;

		// read the header
		size_t headread = fread( &header, sizeof(cyHairFileHeader), 1, fp );

		#define _CY_FAILED_RETURN(errorno) { Initialize(); fclose( fp ); return errorno; }


		// Check if header is correctly read
		if ( headread < 1 ) _CY_FAILED_RETURN(CY_HAIR_FILE_ERROR_CANT_READ_HEADER);

		// Check if this is a hair file
		if ( strncmp( header.signature, "HAIR", 4) != 0 ) _CY_FAILED_RETURN(CY_HAIR_FILE_ERROR_WRONG_SIGNATURE);

		// Read helper
		auto read = [&]<typename T>(
			std::vector<T>& arr, unsigned int mask, size_t dimension
		) {
			if ( !(header.arrays & mask)) return 0;
			arr.clear();
			// Segments are stored per-hair; all other arrays are per-point
			const size_t count = (mask == CY_HAIR_FILE_SEGMENTS_BIT) ?
				header.hair_count : header.point_count * dimension;
			arr.resize(count);
			size_t readcount = fread(
				arr.data(), sizeof(T), count, fp
			);
			if ( readcount < count )
				_CY_FAILED_RETURN(CY_HAIR_FILE_ERROR_READING_SEGMENTS);
			return 0;
		};

		// Read arrays
		read(segments, CY_HAIR_FILE_SEGMENTS_BIT, 1);
		read(points, CY_HAIR_FILE_POINTS_BIT, 3);
		read(thickness, CY_HAIR_FILE_THICKNESS_BIT, 1);
		read(transparency, CY_HAIR_FILE_TRANSPARENCY_BIT, 1);
		read(colors, CY_HAIR_FILE_COLORS_BIT, 3);
		read(uvs, CY_HAIR_FILE_UVS_BIT, 2);

		fclose( fp );

		return header.hair_count;
	}

	/// Saves hair data to the given HAIR file.
	int SaveToFile( const char *filename ) const
	{
		FILE *fp;
		fp = fopen( filename, "wb" );
		if ( fp == NULL ) return -1;

		// Write header
		fwrite( &header, sizeof(cyHairFileHeader), 1, fp );

		// Write arrays
		if ( header.arrays & CY_HAIR_FILE_SEGMENTS_BIT ) fwrite( segments.data(), sizeof(unsigned short), header.hair_count, fp );
		if ( header.arrays & CY_HAIR_FILE_POINTS_BIT ) fwrite( points.data(), sizeof(float), header.point_count*3, fp );
		if ( header.arrays & CY_HAIR_FILE_THICKNESS_BIT ) fwrite( thickness.data(), sizeof(float), header.point_count, fp );
		if ( header.arrays & CY_HAIR_FILE_TRANSPARENCY_BIT ) fwrite( transparency.data(), sizeof(float), header.point_count, fp );
		if ( header.arrays & CY_HAIR_FILE_COLORS_BIT ) fwrite( colors.data(), sizeof(float), header.point_count*3, fp );
		if ( header.arrays & CY_HAIR_FILE_UVS_BIT ) fwrite( uvs.data(), sizeof(float), header.point_count*2, fp );

		fclose( fp );

		return header.hair_count;
	}


	//////////////////////////////////////////////////////////////////////////
	///@name Other Methods

	/// Fills the given direction array with normalized directions using the points array.
	/// Call this function if you need strand directions for shading.
	/// The given array dir should be allocated as an array of size 3 times point count.
	/// Returns point count, returns zero if fails.
	int FillDirectionArray( float *dir )
	{
		if ( dir==NULL || header.point_count<=0 || points.empty() ) return 0;

		int p = 0;	// point index
		for ( unsigned int i=0; i<header.hair_count; i++ ) {
			int s = (not segments.empty()) ? segments[i] : header.d_segments;
			if ( s > 1 ) {
				// direction at point1
				float len0, len1;
				ComputeDirection(
					span3f(&dir[(p + 1) * 3], 3),
					len0,
					len1,
					span3f(&points[p * 3], 3),
					span3f(&points[(p + 1) * 3], 3),
					span3f(&points[(p + 2) * 3], 3)
				);

				// direction at point0
				float d0[3];
				d0[0] = points[(p+1)*3]   - dir[(p+1)*3]  *len0*0.3333f - points[p*3];
				d0[1] = points[(p+1)*3+1] - dir[(p+1)*3+1]*len0*0.3333f - points[p*3+1];
				d0[2] = points[(p+1)*3+2] - dir[(p+1)*3+2]*len0*0.3333f - points[p*3+2];
				float d0lensq = d0[0]*d0[0] + d0[1]*d0[1] + d0[2]*d0[2];
				float d0len = ( d0lensq > 0 ) ? (float) sqrt(d0lensq) : 1.0f;
				dir[p*3]   = d0[0] / d0len;
				dir[p*3+1] = d0[1] / d0len;
				dir[p*3+2] = d0[2] / d0len;

				// We computed the first 2 points
				p += 2;

				// Compute the direction for the rest
				for ( int t=2; t<s; t++, p++ ) {
					ComputeDirection(
						span3f(&dir[p * 3], 3),
						len0, len1,
						span3f(&points[(p - 1) * 3], 3),
						span3f(&points[p * 3], 3),
						span3f(&points[(p + 1) * 3], 3)
					);
				}

				// direction at the last point
				d0[0] = - points[(p-1)*3]   + dir[(p-1)*3]  *len1*0.3333f + points[p*3];
				d0[1] = - points[(p-1)*3+1] + dir[(p-1)*3+1]*len1*0.3333f + points[p*3+1];
				d0[2] = - points[(p-1)*3+2] + dir[(p-1)*3+2]*len1*0.3333f + points[p*3+2];
				d0lensq = d0[0]*d0[0] + d0[1]*d0[1] + d0[2]*d0[2];
				d0len = ( d0lensq > 0 ) ? (float) sqrt(d0lensq) : 1.0f;
				dir[p*3]   = d0[0] / d0len;
				dir[p*3+1] = d0[1] / d0len;
				dir[p*3+2] = d0[2] / d0len;
				p++;

			} else if ( s > 0 ) {
				// if it has a single segment
				float d0[3];
				d0[0] = points[(p+1)*3]   - points[p*3];
				d0[1] = points[(p+1)*3+1] - points[p*3+1];
				d0[2] = points[(p+1)*3+2] - points[p*3+2];
				float d0lensq = d0[0]*d0[0] + d0[1]*d0[1] + d0[2]*d0[2];
				float d0len = ( d0lensq > 0 ) ? (float) sqrt(d0lensq) : 1.0f;
				dir[p*3]   = d0[0] / d0len;
				dir[p*3+1] = d0[1] / d0len;
				dir[p*3+2] = d0[2] / d0len;
				dir[(p+1)*3]   = dir[p*3];
				dir[(p+1)*3+1] = dir[p*3+1];
				dir[(p+1)*3+2] = dir[p*3+2];
				p += 2;
			}
			//*/
		}
		return p;
	}


private:
	//////////////////////////////////////////////////////////////////////////
	///@name Private Variables and Methods

	cyHairFileHeader header;
	std::vector<unsigned short>	segments;
	std::vector<float>			points;
	std::vector<float>			thickness;
	std::vector<float>			transparency;
	std::vector<float>			colors;
	std::vector<float>			uvs;

	using span3f = std::span<float, 3>;

	// Given point before (p0) and after (p2), computes the direction (d) at p1.
	float ComputeDirection(span3f d, float &d0len, float &d1len, const span3f p0, const span3f p1, const span3f p2 )
	{
		// line from p0 to p1
		float d0[3];
		d0[0] = p1[0] - p0[0];
		d0[1] = p1[1] - p0[1];
		d0[2] = p1[2] - p0[2];
		float d0lensq = d0[0]*d0[0] + d0[1]*d0[1] + d0[2]*d0[2];
		d0len = ( d0lensq > 0 ) ? (float) sqrt(d0lensq) : 1.0f;

		// line from p1 to p2
		float d1[3];
		d1[0] = p2[0] - p1[0];
		d1[1] = p2[1] - p1[1];
		d1[2] = p2[2] - p1[2];
		float d1lensq = d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2];
		d1len = ( d1lensq > 0 ) ? (float) sqrt(d1lensq) : 1.0f;

		// make sure that d0 and d1 has the same length
		d0[0] *= d1len / d0len;
		d0[1] *= d1len / d0len;
		d0[2] *= d1len / d0len;

		// direction at p1
		d[0] = d0[0] + d1[0];
		d[1] = d0[1] + d1[1];
		d[2] = d0[2] + d1[2];
		float dlensq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
		float dlen = ( dlensq > 0 ) ? (float) sqrt(dlensq) : 1.0f;
		d[0] /= dlen;
		d[1] /= dlen;
		d[2] /= dlen;

		return d0len;
	}
};

//-------------------------------------------------------------------------------

namespace cy {
	typedef cyHairFileHeader HairFileHeader;
	typedef cyHairFile HairFile;
}

//-------------------------------------------------------------------------------

}

#endif
// vim: autoindent noexpandtab tabstop=4 shiftwidth=4
