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
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/GeometryVisitor.h>
#include <SFCGAL/algorithm/normal.h>

namespace SFCGAL {

///
///
///
PolyhedralSurface::PolyhedralSurface():
	Surface(),
	_polygons()
{

}

///
///
///
PolyhedralSurface::PolyhedralSurface( const std::vector< Polygon > & polygons ) :
	Surface()
{
	for ( size_t i = 0; i < polygons.size(); i++ ){
		_polygons.push_back( polygons[i].clone() ) ;
	}
}

///
///
///
PolyhedralSurface::PolyhedralSurface( PolyhedralSurface const& other ) :
	Surface(),
	_polygons(other._polygons)
{

}

//
// Segment ring type
// Used for simplification
class Ring
{
public:
	typedef CGAL::Segment_3<Kernel> value_type;
	typedef std::list<value_type> Segments;
	typedef Segments::const_iterator const_iterator;
	typedef Segments::iterator iterator;
	typedef CGAL::Circulator_from_container<Segments> Circulator;
	typedef CGAL::Const_circulator_from_container<Segments> Const_circulator;

	Ring() : segments(), insert_location( segments.end() ) {
	}

	Ring( const Ring& other ) {
		// insert_location must point to the equivalent copied list element
		segments = other.segments;
		int idx = std::distance( other.segments.begin(), (const_iterator)other.insert_location );
		insert_location = segments.begin();
		std::advance( insert_location, idx );
	}

	iterator begin() { return segments.begin(); }
	iterator end() { return segments.end(); }
	const_iterator begin() const { return segments.begin(); }
	const_iterator end() const { return segments.end(); }

	// returns a circulator starting at it
	Circulator circulator( iterator it ) { return Circulator( &segments, it ); }
	Const_circulator const_circulator( const_iterator it ) const { return Const_circulator( &segments, it ); }

	void push( const value_type& v ) {
		Segments::iterator it = segments.insert( insert_location, v );
		insert_location = ++it;
	}

	iterator erase( iterator it ) {
		insert_location = segments.erase( it );
		return insert_location;
	}

	iterator back() {
		return insert_location;
	}

	void clear() {
		segments.clear();
		insert_location = segments.end();
	}

	/**
	 * Merges a given ring
	 * Returns true if the merge was possible
	 */
	bool merge( const Ring& other ) {
		// common edge on this
		Segments::iterator common_edge = segments.end();
		// common edge on other
		Segments::const_iterator common_edge2 = other.end();

		if ( segments.empty() ) {
			for ( Segments::const_iterator it = other.begin(); it != other.end(); ++it ) {
				this->push( *it );
			}
			return true;
		}

		for ( Segments::const_iterator it = other.begin(); it != other.end(); ++it ) {
			common_edge = std::find( segments.begin(), segments.end(), it->opposite() );
			if ( common_edge != segments.end() ) {
				this->erase( common_edge );
				common_edge2 = it;
				break;
			}
		}
		if ( common_edge2 == other.end() ) {
			return false;
		}

		// circulator around facet, starting after common_edge2
		Const_circulator cit = other.const_circulator( common_edge2 );
		Const_circulator cit_end = other.const_circulator( common_edge2 );
		++cit;
		do {
				this->push( *cit );
				++cit;
		} while ( cit != cit_end );

		return true;
	}

	/**
	 * Convert a ring to a Polygon
	 * Filter out collinear points
	 */
	SFCGAL::Polygon toPolygon() const
	{
		LineString boundary;

		for ( Ring::const_iterator rit = segments.begin(); rit != segments.end(); ++rit ) {
			// filter out collinear points
			Ring::const_iterator next;
			next = rit;
			next++;
			if ( next == segments.end() ) {
				next = segments.begin();
			}
			
			if ( CGAL::collinear( rit->source(), rit->target(), next->target() ) ) {
				continue;
			}
			boundary.addPoint( rit->target() );
		}
		boundary.addPoint( boundary.pointN(0) );
		
		return Polygon( boundary );
	}

private:
	Segments segments;
	// where to insert the next segment ?
	Segments::iterator insert_location;
};

struct SurfaceProcessing
{
	typedef std::list<Ring> Rings;
	Ring merged;
	Rings not_merged;
};

///
///
///
PolyhedralSurface::PolyhedralSurface( const MarkedPolyhedron& poly, bool simplify ) :
Surface()
{
	if ( ! simplify ) {

		for ( MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin(); fit != poly.facets_end(); ++fit ) {
			LineString* face = new LineString();
			MarkedPolyhedron::Halfedge_around_facet_const_circulator hit = fit->facet_begin();
			do {
				face->addPoint( hit->vertex()->point() );
				++hit;
			} while ( hit != fit->facet_begin() );
			// close the ring
			face->addPoint( hit->vertex()->point() );
			_polygons.push_back( new Polygon( face ) );
		}
		return;
	}
	// else
	
	typedef CGAL::Plane_3<Kernel> Plane;
	typedef std::list<std::pair<Plane, SurfaceProcessing> > Plane2Surfaces;
	
	Plane2Surfaces plane2Surfaces;
	
	// compute planes
	for ( MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin(); fit != poly.facets_end(); ++fit ) {
		std::vector<CGAL::Point_3<Kernel> > pts;
		Ring ring;
		MarkedPolyhedron::Halfedge_around_facet_const_circulator hit = fit->facet_begin();
		do {
			pts.push_back( hit->vertex()->point() );
			const CGAL::Point_3<Kernel>& s = hit->opposite()->vertex()->point();
			const CGAL::Point_3<Kernel>& t = hit->vertex()->point();
			ring.push( CGAL::Segment_3<Kernel>( s, t ) );
			hit++;
		} while ( hit != fit->facet_begin() );
		pts.push_back( hit->vertex()->point() );
		
		CGAL::Vector_3<Kernel> vn = SFCGAL::algorithm::normal3D( pts );
		CGAL::Plane_3<Kernel> plane( fit->facet_begin()->vertex()->point(), vn );
		
		// insert the new plane if needed
		Plane2Surfaces::iterator pit = plane2Surfaces.end();
		for ( pit = plane2Surfaces.begin(); pit != plane2Surfaces.end(); ++pit ) {
			// use the Plane::operator==
			if ( pit->first == plane ) {
				break;
			}
		}
		if ( pit == plane2Surfaces.end() ) {
			// no plane found, have to insert it
			pit = plane2Surfaces.insert( plane2Surfaces.end(), std::make_pair( plane, SurfaceProcessing() ) );
		}
		
		SurfaceProcessing& surface = pit->second;
		// this facet is to be processed
		surface.not_merged.push_back( ring );
	}

	// for each plane
	for ( Plane2Surfaces::iterator it = plane2Surfaces.begin(); it != plane2Surfaces.end(); ++it ) {
		// for each facet to be processed
		SurfaceProcessing& surface = it->second;
		while ( ! surface.not_merged.empty() ) {
			bool merged = false;
			for ( SurfaceProcessing::Rings::iterator rit = surface.not_merged.begin();
			      rit != surface.not_merged.end();
			      ++rit ) {
				// try to merge the facet into the ring
				if ( surface.merged.merge( *rit ) ) {
					merged = true;
					surface.not_merged.erase( rit );
					break;
				}
			}
			if ( !merged || surface.not_merged.empty() ) {
				addPolygon( surface.merged.toPolygon() );
				surface.merged.clear();
			}
		}
	}
}

///
///
///
PolyhedralSurface& PolyhedralSurface::operator = ( const PolyhedralSurface & other )
{
	_polygons = other._polygons ;
	return *this ;
}


///
///
///
PolyhedralSurface::~PolyhedralSurface()
{

}

///
///
///
PolyhedralSurface * PolyhedralSurface::clone() const
{
	return new PolyhedralSurface(*this);
}

///
///
///
std::string  PolyhedralSurface::geometryType() const
{
	return "PolyhedralSurface" ;
}

///
///
///
GeometryType  PolyhedralSurface::geometryTypeId() const
{
	return TYPE_POLYHEDRALSURFACE ;
}

///
///
///
int  PolyhedralSurface::dimension() const
{
	return 2 ;
}

///
///
///
int  PolyhedralSurface::coordinateDimension() const
{
	if ( isEmpty() ){
		return 0 ;
	}else{
		return _polygons.front().coordinateDimension() ;
	}
}

///
///
///
bool  PolyhedralSurface::isEmpty() const
{
	return _polygons.empty();
}

///
///
///
bool  PolyhedralSurface::is3D() const
{
	if ( isEmpty() ){
		return false ;
	}else{
		return _polygons.front().is3D() ;
	}
}

///
///
///
bool  PolyhedralSurface::isMeasured() const
{
	if ( isEmpty() ){
		return false ;
	}else{
		return _polygons.front().isMeasured() ;
	}
}



///
///
///
TriangulatedSurface  PolyhedralSurface::toTriangulatedSurface() const
{
	TriangulatedSurface result ;
	triangulate::triangulatePolygon3D( *this, result );
	return result ;
}

///
///
///
void  PolyhedralSurface::addPolygon( const Polygon & polygon )
{
	addPolygon( polygon.clone() );
}

///
///
///
void  PolyhedralSurface::addPolygon( Polygon* polygon )
{
	BOOST_ASSERT( polygon != NULL );
	_polygons.push_back( polygon );
}

///
///
///
void  PolyhedralSurface::addPolygons( const PolyhedralSurface & polyhedralSurface )
{
	for ( size_t i = 0; i < polyhedralSurface.numPolygons(); i++ ){
		addPolygon( polyhedralSurface.polygonN(i) );
	}
}

///
///
///
size_t   PolyhedralSurface::numGeometries() const
{
	return _polygons.size() ;
}


///
///
///
const Polygon  &   PolyhedralSurface::geometryN( size_t const& n ) const
{
	return _polygons[n] ;
}

///
///
///
Polygon & PolyhedralSurface::geometryN( size_t const& n )
{
	return _polygons[n];
}

///
///
///
void PolyhedralSurface::accept( GeometryVisitor & visitor )
{
	return visitor.visit(*this);
}

///
///
///
void PolyhedralSurface::accept( ConstGeometryVisitor & visitor ) const
{
	return visitor.visit(*this);
}
}//SFCGAL



