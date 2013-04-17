/**
 *   SFCGAL
 *
 *   Copyright (C) 2012-2013 Oslandia <infos@oslandia.com>
 *   Copyright (C) 2012-2013 IGN (http://www.ign.fr)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <SFCGAL/algorithm/difference.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/tools/Registry.h>
#include <SFCGAL/algorithm/area.h>
#include <SFCGAL/Envelope.h>
#include <SFCGAL/io/GeometryStreams.h>
#include <SFCGAL/all.h>
#include <SFCGAL/transform/AffineTransform3.h>

#include <SFCGAL/io/wkt.h>

#include <boost/test/unit_test.hpp>

using namespace SFCGAL;
using namespace boost::unit_test ;

BOOST_AUTO_TEST_SUITE( SFCGAL_algorithm_DifferenceTest )

#if 0
BOOST_AUTO_TEST_CASE( testDifferenceXPoint )
{
	// The same point
	BOOST_CHECK( algorithm::difference( Point(0,0), Point(0,0) )->isEmpty() );
	// A different point
	BOOST_CHECK( *algorithm::difference( Point(1,0), Point(0,0) ) == Point(1,0) );

	// check difference(X, point) == X
	std::vector<std::string> typeNames;
	for ( size_t i = 0; i < typeNames.size(); ++i ) {
		std::auto_ptr<Geometry> newGeo( tools::Registry::instance().newGeometryByTypeName( typeNames[i] ) );
		std::auto_ptr<Geometry> diffGeo = algorithm::difference( *newGeo, Point(0, 0) );
		BOOST_CHECK( *newGeo == *diffGeo );
	}
}

BOOST_AUTO_TEST_CASE( testDifference3DXPoint )
{
	// The same point
	BOOST_CHECK( algorithm::difference3D( Point(0,0,0), Point(0,0,0) )->isEmpty() );
	// A different point
	BOOST_CHECK( *algorithm::difference3D( Point(1,0,0), Point(0,0,0) ) == Point(1,0,0) );

	// check difference(X, point) == X
	std::vector<std::string> typeNames;
	for ( size_t i = 0; i < typeNames.size(); ++i ) {
		std::auto_ptr<Geometry> newGeo( tools::Registry::instance().newGeometryByTypeName( typeNames[i] ) );
		std::auto_ptr<Geometry> diffGeo = algorithm::difference3D( *newGeo, Point(0, 0, 0) );
		BOOST_CHECK( *newGeo == *diffGeo );
	}
}

BOOST_AUTO_TEST_CASE( testDifferenceXLineString )
{
	// Point x Linestring intersecting
	BOOST_CHECK( algorithm::difference( Point(0,0), *io::readWkt("LINESTRING(0 0,1 1)"))->isEmpty() );
	// Point x linestring not intersecting
	BOOST_CHECK( *algorithm::difference( Point(0,0), *io::readWkt("LINESTRING(0 1,1 1)")) == Point(0, 0) );

	// two linestrings with two segments overlapping
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0.5 0,0.7 0)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("MULTILINESTRING((0 0,0.5 0),(0.7 0,1 0))") );
	}
	// two linestrings with two opposite segments overlapping
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0.7 0,0.5 0)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("MULTILINESTRING((0 0,0.5 0),(0.7 0,1 0))") );
	}
	// two linestrings with two segments crossing
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(-1 0,1 0)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0 -1,0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( *diff == *ls1 );
	}
	// two linestrings with two segments partly overlapping
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(-1 0,0.7 0)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(0.7 0,1 0)") );
	}
	// two linestrings with a segment covering another one
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(-1 0,2 0)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( diff->isEmpty() );
	}
	// two linestrings that do not intersect
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0 1,1 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( *diff == *ls1 );
	}
	// two linestrings with more than one segment
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0,1 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0.3 0,1 0,1 0.4)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("MULTILINESTRING((0 0,0.3 0),(1 0.4,1 1))") );
	}

	// a linestring and a surface
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0,1 0)");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((0.5 0,2 0,2 1,0.5 1,0.5 0))");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls1, *square2 );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(0 0,0.5 0)") );
	}

	// check difference(X, linestring) == X, with dimension(X) > 1
	// TODO: add generators of random geometries to avoid empty geometries here ?
	std::vector<std::string> typeNames;
	for ( size_t i = 0; i < typeNames.size(); ++i ) {
		std::auto_ptr<Geometry> newGeo( tools::Registry::instance().newGeometryByTypeName( typeNames[i] ) );
		if ( newGeo->dimension() > 1 ) {
			std::auto_ptr<Geometry> diffGeo = algorithm::difference( *newGeo, Point(0, 0) );
			BOOST_CHECK( *newGeo == *diffGeo );
		}
	}
}

BOOST_AUTO_TEST_CASE( testDifference3DXLineString )
{
	// Point x Linestring intersecting
	BOOST_CHECK( algorithm::difference3D( Point(0,0,0), *io::readWkt("LINESTRING(0 0 0,1 1 0)"))->isEmpty() );
	// Point x linestring not intersecting
	BOOST_CHECK( *algorithm::difference3D( Point(0,0,0), *io::readWkt("LINESTRING(0 1 0,1 1 0)")) == Point(0, 0, 0) );

	// two linestrings with two segments overlapping
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0.5 0 1,0.7 0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("MULTILINESTRING((0 0 1,0.5 0 1),(0.7 0 1,1 0 1))") );
	}
	// two linestrings with two opposite segments overlapping
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0.7 0 1,0.5 0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("MULTILINESTRING((0 0 1,0.5 0 1),(0.7 0 1,1 0 1))") );
	}
	// two linestrings with two segments crossing
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(-1 0 1,1 0 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0 -1 1,0 1 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( *diff == *ls1 );
	}
	// two linestrings with two segments partly overlapping
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(-1 0 1,0.7 0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(0.7 0 1,1 0 1)") );
	}
	// two linestrings with a segment covering another one
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(-1 0 1,2 0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( diff->isEmpty() );
	}
	// two linestrings that do not intersect
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0 1 1,1 1 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( *diff == *ls1 );
	}
	// two linestrings with more than one segment
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1,1 1 1)");
		std::auto_ptr<Geometry> ls2 = io::readWkt("LINESTRING(0.3 0 1,1 0 1,1 0.4 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *ls2 );
		BOOST_CHECK( *diff == *io::readWkt("MULTILINESTRING((0 0 1,0.3 0 1),(1 0.4 1,1 1 1))") );
	}

	// a linestring and a surface
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,1 0 1)");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((0.5 0 1,2 0 1,2 1 1,0.5 1 1,0.5 0 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *square2 );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(0 0 1,0.5 0 1)") );
	}

	// a linestring and a volume
	{
		std::auto_ptr<Geometry> ls1 = io::readWkt("LINESTRING(0 0 1,2 0 1)");
		std::auto_ptr<Geometry> cube2 = io::readWkt("SOLID(("
							    "((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0))," // front face
							    "((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0))," // right face
							    "((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))," // top face
							    "((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1))," // left face
							    "((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1))," // back face
							    "((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))" // bottom face
							    "))");

		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls1, *cube2 );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(1 0 1,2 0 1)") );
	}

	// check difference(X, linestring) == X, with dimension(X) > 1
	// TODO: add generators of random geometries to avoid empty geometries here ?
	std::vector<std::string> typeNames;
	for ( size_t i = 0; i < typeNames.size(); ++i ) {
		std::auto_ptr<Geometry> newGeo( tools::Registry::instance().newGeometryByTypeName( typeNames[i] ) );
		if ( newGeo->dimension() > 1 ) {
			std::auto_ptr<Geometry> diffGeo = algorithm::difference3D( *newGeo, Point(0, 0) );
			BOOST_CHECK( *newGeo == *diffGeo );
		}
	}
}

BOOST_AUTO_TEST_CASE( testDifferenceXPolygon )
{
	// Point x polygon intersecting
	BOOST_CHECK( algorithm::difference( Point(0,0), *io::readWkt("POLYGON((0 0,1 0,1 1,0 1,0 0))"))->isEmpty() );
	// Point x polygon not intersecting
	BOOST_CHECK( *algorithm::difference( Point(-1,0), *io::readWkt("POLYGON((0 0,1 0,1 1,0 1,0 0))")) == Point(-1, 0) );

	{
		// a linestring partly cut by a square
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((0 0,1 0,1 1,0 1,0 0))");
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(-1 0,1 0)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls, *square );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(-1 0,0 0)") );
	}
	{
		// a linestring outside a square
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((0 0,1 0,1 1,0 1,0 0))");
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(-2 0,-1 0)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls, *square );
		BOOST_CHECK( *diff == *ls );
	}
	{
		// a linestring partly inside the hole of a square
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0),(2 2,2 8,8 8,8 2,2 2))");
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(1 1,3 3)");
		std::auto_ptr<Geometry> diff = algorithm::difference( *ls, *square );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(2 2,3 3)") );
	}

	{
		// a polygon touching another polygon
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((10 0,20 0,20 10,10 10,10 0))");
		std::auto_ptr<Geometry> diff = algorithm::difference( *square1, *square2 );
		BOOST_CHECK( *diff == *square1 );
	}
	{
		// a polygon partly crossing another polygon
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((5 0,15 0,15 10,5 10,5 0))");
		std::auto_ptr<Geometry> diff = algorithm::difference( *square1, *square2 );
		BOOST_CHECK( *diff == *io::readWkt("POLYGON((0 0,5 0,5 10,0 10,0 0))") );
	}
	{
		// a polygon partly crossing another polygon in a hole
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0),(2 2,2 8,8 8,8 2,2 2))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((1 1,5 1,5 5,1 5,1 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference( *square1, *square2 );
		BOOST_CHECK( *diff ==  *io::readWkt("POLYGON((0 0,10 0,10 10,0 10,0 0),(8 8,8 2,5 2,5 1,1 1,1 5,2 5,2 8,8 8))") );
	}
	{
		// two polygons that do not intersect
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((5 0,10 0,10 10,5 10,5 0))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((0 0,2 0,2 2,0 2,0 0))");
		std::auto_ptr<Geometry> diff = algorithm::difference( *square1, *square2 );
		BOOST_CHECK( *diff == *square1 );
	}
}

BOOST_AUTO_TEST_CASE( testDifference3DXPolygon )
{
	// Point x polygon intersecting
	BOOST_CHECK( algorithm::difference3D( Point(0,0,0), *io::readWkt("POLYGON((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0))"))->isEmpty() );
	// Point x polygon not intersecting
	BOOST_CHECK( *algorithm::difference3D( Point(-1,0,0), *io::readWkt("POLYGON((0 0 0,1 0 0,1 1 0,0 1 0,0 0 0))")) == Point(-1, 0, 0) );

	{
		// a linestring partly cut by a square
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1))");
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(-1 0 1,1 0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls, *square );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(-1 0 1,0 0 1)") );
	}
	{
		// a linestring outside a square
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1))");
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(-2 0 1,-1 0 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls, *square );
		BOOST_CHECK( *diff == *ls );
	}
	{
		// a linestring partly inside the hole of a square
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((0 0 1,10 0 1,10 10 1,0 10 1,0 0 1),(2 2 1,2 8 1,8 8 1,8 2 1,2 2 1))");
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(1 1 1,3 3 1)");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *ls, *square );
		BOOST_CHECK( *diff == *io::readWkt("LINESTRING(2 2 1,3 3 1)") );
	}

	{
		// a triangle crossing a polygon
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0 1,10 0 1,10 10 1,0 10 1,0 0 1))");
		std::auto_ptr<Geometry> tri2 = io::readWkt("TRIANGLE((0 0 1,10 0 1,10 10 1,0 0 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square1, *tri2 );
		BOOST_CHECK( diff->envelope() == Envelope( 0, 10, 0, 10, 1, 1 ) );
		BOOST_CHECK( algorithm::area3D( *diff ) == 50 );
	}
	{
		// a triangle crossing a polygon
		std::auto_ptr<Geometry> tri1 = io::readWkt("TRIANGLE((0 0 1,10 0 1,0 10 1,0 0 1))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((5 0 1,15 0 1,15 10 1,5 10 1,5 0 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *tri1, *square2 );
		BOOST_CHECK( diff->envelope() == Envelope( 0, 5, 0, 10, 1, 1 ) );
		BOOST_CHECK( algorithm::area3D( *diff ) == 37.5 );
	}

	{
		// a polygon touching another polygon
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0 1,10 0 1,10 10 1,0 10 1,0 0 1))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((10 0 1,20 0 1,20 10 1,10 10 1,10 0 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square1, *square2 );
		// diff == square1
		BOOST_CHECK( diff->envelope() == Envelope( 0, 10, 0, 10, 1, 1 ) );
		BOOST_CHECK( algorithm::area3D( *diff ) == 100 );
	}
	{
		// a polygon partly crossing another polygon
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0 1,10 0 1,10 10 1,0 10 1,0 0 1))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((5 0 1,15 0 1,15 10 1,5 10 1,5 0 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square1, *square2 );
		BOOST_CHECK( diff->envelope() == Envelope( 0, 5, 0, 10, 1, 1 ) );
		BOOST_CHECK( algorithm::area3D( *diff ) == 50 );
	}
	{
		// a polygon partly crossing another polygon in a hole
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0 1,10 0 1,10 10 1,0 10 1,0 0 1),(2 2 1,2 8 1,8 8 1,8 2 1,2 2 1))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((1 1 1,5 1 1,5 5 1,1 5 1,1 1 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square1, *square2 );
		BOOST_CHECK( diff->envelope() == Envelope( 0, 10, 0, 10, 1, 1 ) );
		BOOST_CHECK( algorithm::area3D( *diff ) == 100 - 6*6 - 4 - 3 );
	}
	{
		// two polygons that do not intersect
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((5 0 1,10 0 1,10 10 1,5 10 1,5 0 1))");
		std::auto_ptr<Geometry> square2 = io::readWkt("POLYGON((0 0 1,2 0 1,2 2 1,0 2 1,0 0 1))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square1, *square2 );
		BOOST_CHECK( diff->envelope() == Envelope( 5, 10, 0, 10, 1, 1 ) );
		BOOST_CHECK( algorithm::area3D( *diff ) == 50 );
	}

	{
		// a surface and a volume
		std::auto_ptr<Geometry> square1 = io::readWkt("POLYGON((0 0 1,10 0 1,10 10 1,0 10 1,0 0 1))");
		std::auto_ptr<Geometry> cube2 = io::readWkt("SOLID(("
							    "((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0))," // front face
							    "((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0))," // right face
							    "((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))," // top face
							    "((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1))," // left face
							    "((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1))," // back face
							    "((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))" // bottom face
							    "))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square1, *cube2 );
		BOOST_CHECK( diff->envelope() == Envelope( 0, 10, 0, 10, 1, 1 ) );
		double a = algorithm::area3D( *diff );
		BOOST_CHECK( (a - 99.0) * (a - 99.0) < 0.001 );
	}
}
#endif

BOOST_AUTO_TEST_CASE( testDifference3DXSolid )
{
	std::auto_ptr<Geometry> cube = io::readWkt("SOLID(("
						   "((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0))," // front face
						   "((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0))," // right face
						   "((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))," // top face
						   "((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1))," // left face
						   "((1 0 1,1 1 1,0 1 1,0 0 1,1 0 1))," // back face
						   "((1 0 0,1 0 1,0 0 1,0 0 0,1 0 0))" // bottom face
						   "))");
	// Point - Solid intersecting
	BOOST_CHECK( algorithm::difference3D( Point(0,0,0), *cube )->isEmpty() );
	// Point - Solid not intersecting
	BOOST_CHECK( *algorithm::difference3D( Point(-1,0,0), *cube) == Point(-1, 0, 0) );

	{
		// linestring - Solid
		std::auto_ptr<Geometry> ls = io::readWkt("LINESTRING(-1 0 0,0.5 0 0)");
		BOOST_CHECK( *algorithm::difference3D( *ls, *cube ) == *io::readWkt("LINESTRING(-1 0 0,0 0 0)") );

		// a linestring partly inside the hole of a cube
		// TODO
	}

	{
		// Triangle - Solid
		std::auto_ptr<Geometry> tri1 = io::readWkt("TRIANGLE((-1 0 0,0.5 0 0,0.5 1 0,-1 0 0))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *tri1, *cube );
		BOOST_CHECK( *diff == *io::readWkt("TRIANGLE((0/1 2/3 0/1,-1/1 0/1 0/1,0/1 0/1 0/1,0/1 2/3 0/1))") );
	}

	{
		// Polygon - Solid
		std::auto_ptr<Geometry> square = io::readWkt("POLYGON((-1 0 0,0.5 0 0,0.5 1 0,-1 1 0,-1 0 0))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *square, *cube );
		double ae = algorithm::area3D( *diff ) - 1.0;
		BOOST_CHECK( ae*ae < 0.001 );
		BOOST_CHECK( diff->envelope() == Envelope( -1, 0, 0, 1, 0, 0 ) );
	}

	{
		// Solid - Solid

		// translation of the cube
		Solid cube1( cube->as<Solid>() );
		CGAL::Vector_3<Kernel> tv( CGAL::Point_3<Kernel>( 0.0, 0.0, 0.0 ),
					   CGAL::Point_3<Kernel>( 0.5, 0.0, 0.0 ));
		transform::AffineTransform3<Kernel> t( CGAL::Aff_transformation_3<Kernel>( CGAL::Translation(),
											   tv ));
		t.transform( cube1 );
		
		std::auto_ptr<Geometry> diff = algorithm::difference3D( cube1, *cube );
		BOOST_CHECK( *diff == *io::readWkt("SOLID((((1 1 0,1 1 1,1.5 1 1,1.5 1 0,1 1 0)),((1.5 1 1,1.5 0 1,1.5 0 0,1.5 1 0,1.5 1 1)),((1 1 0,1 0 0,1 0 1,1 1 1,1 1 0)),((1.5 1 0,1.5 0 0,1 0 0,1 1 0,1.5 1 0)),((1.5 1 1,1 1 1,1 0 1,1.5 0 1,1.5 1 1)),((1.5 0 1,1 0 1,1 0 0,1.5 0 0,1.5 0 1))))") );
	}

	{
		// Solid - Solid

		std::auto_ptr<Geometry> cube2 = io::readWkt("SOLID(("
							    "((0.2 0 0,0.2 1 0,0.8 1 0,0.8 0 0,0.2 0 0))," // front face
							    "((0.8 0 0,0.8 1 0,0.8 1 1,0.8 0 1,0.8 0 0))," // right face
							    "((0.2 1 0,0.2 1 1,0.8 1 1,0.8 1 0,0.2 1 0))," // top face
							    "((0.2 0 1,0.2 1 1,0.2 1 0,0.2 0 0,0.2 0 1))," // left face
							    "((0.8 0 1,0.8 1 1,0.2 1 1,0.2 0 1,0.8 0 1))," // back face
							    "((0.8 0 0,0.8 0 1,0.2 0 1,0.2 0 0,0.8 0 0))" // bottom face
							    "))");
		std::auto_ptr<Geometry> diff = algorithm::difference3D( *cube, *cube2 );
		BOOST_CHECK( *diff == *io::readWkt("SOLID((((0.2 0 1,0 0 1,0 0 0,0.2 0 0,0.2 0 1)),((1 0 1,0.8 0 1,0.8 0 0,1 0 0,1 0 1)),((0.2 0 0,0.2 1 0,0.2 1 1,0.2 0 1,0.2 0 0)),((0 0 0,0 1 0,0.2 1 0,0.2 0 0,0 0 0)),((0.8 0 0,0.8 1 0,1 1 0,1 0 0,0.8 0 0)),((0 0 1,0 1 1,0 1 0,0 0 0,0 0 1)),((0.2 1 1,0 1 1,0 0 1,0.2 0 1,0.2 1 1)),((1 0 1,1 1 1,0.8 1 1,0.8 0 1,1 0 1)),((0.2 1 1,0.2 1 0,0 1 0,0 1 1,0.2 1 1)),((1 1 0,0.8 1 0,0.8 1 1,1 1 1,1 1 0)),((0.8 0 1,0.8 1 1,0.8 1 0,0.8 0 0,0.8 0 1)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0))))") );
	}
}

BOOST_AUTO_TEST_SUITE_END()

