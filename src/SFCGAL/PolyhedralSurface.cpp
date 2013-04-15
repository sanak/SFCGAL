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

	Ring() : segments(), insert_location( segments.end() ) {
	}

	Ring( const Ring& other ) {
		segments = other.segments;
		int idx = std::distance( other.segments.begin(), (const_iterator)other.insert_location );
		insert_location = segments.begin();
		for ( int i = 0; i < idx; ++i ) {
			insert_location++;
		}
	}

	iterator begin() { return segments.begin(); }
	iterator end() { return segments.end(); }
	const_iterator begin() const { return segments.begin(); }
	const_iterator end() const { return segments.end(); }

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

private:
	Segments segments;
	// where to insert the next segment ?
	Segments::iterator insert_location;
};

struct FacetAnnotationElement
{
	CGAL::Plane_3<Kernel> plane;
	bool merged;
	FacetAnnotationElement() : merged(false) {}
	FacetAnnotationElement( const CGAL::Plane_3<Kernel>& p ) : plane(p), merged(false) {}
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
	}
	else {
		//
		// Idea of the algorithm:
		// Repeat until all facets are processed:
		//   For each facet Fi:
		//     Compute the plane equation P
		//     Retrieve the ring R associated with the plane P
		//     If it exists:
		//       If Fi does not share an edge Ec with R, continue
		//       Remove Ec from R
		//     Add each other vertices of Fi to R
		//     Mark Fi as processed
		//
		// For each segments of the ring:
		//   add it to the linestring if points are not collinear

		typedef CGAL::Plane_3<Kernel> Plane;
		typedef std::list<std::pair<Plane, Ring> > RingPartMap;

		typedef std::map<MarkedPolyhedron::Facet_const_iterator, FacetAnnotationElement> FacetAnnotation;
		FacetAnnotation facet_annotation;

		RingPartMap ring_parts;

		// for each facet
		bool all_facets_merged = false;

		while ( !all_facets_merged ) {
			all_facets_merged = true;

			for ( MarkedPolyhedron::Facet_const_iterator fit = poly.facets_begin(); fit != poly.facets_end(); ++fit ) {
				FacetAnnotation::iterator ait = facet_annotation.find( fit );
				if ( ait == facet_annotation.end() ) {
					std::vector<CGAL::Point_3<Kernel> > pts;
					MarkedPolyhedron::Halfedge_around_facet_const_circulator hit = fit->facet_begin();
					do {
						pts.push_back( hit->vertex()->point() );
						hit++;
					} while ( hit != fit->facet_begin() );
					pts.push_back( hit->vertex()->point() );

					CGAL::Vector_3<Kernel> vn = SFCGAL::algorithm::normal3D( pts );
					CGAL::Plane_3<Kernel> plane( fit->facet_begin()->vertex()->point(), vn );
					bool ok;
					boost::tie(ait,ok) = facet_annotation.insert( std::make_pair( fit, plane ) );
				}
				// this facet has already been processed, skip
				if ( ait->second.merged ) {
					continue;
				}
				CGAL::Plane_3<Kernel>& plane = ait->second.plane;
				all_facets_merged = false;

				bool new_plane = false;
				RingPartMap::iterator parts;
				for ( parts = ring_parts.begin(); parts != ring_parts.end(); ++parts ) {
					// use Plane::operator==
					if ( parts->first == plane ) {
						break;
					}
				}
				if ( parts == ring_parts.end() ) {
					new_plane = true;
					Ring p;
					parts = ring_parts.insert( ring_parts.end(), std::make_pair( plane, p ) );
				}
				Ring& ring = parts->second;
				
				MarkedPolyhedron::Halfedge_around_facet_const_circulator hit = fit->facet_begin();
				MarkedPolyhedron::Halfedge_around_facet_const_circulator hit_end = fit->facet_begin();
				if ( ! new_plane ) {
					// first pass over vertices : find the shared edge
					Ring::iterator common_edge = ring.end();
					do {
						const CGAL::Point_3<Kernel>& s = hit->opposite()->vertex()->point();
						const CGAL::Point_3<Kernel>& t = hit->vertex()->point();
						
						CGAL::Segment_3<Kernel> opposite_segment( t, s );
						
						// look for the common edge
						common_edge = std::find( ring.begin(), ring.end(), opposite_segment );
						if ( common_edge != ring.end() ) {
							break;
						}
						++hit;
					} while ( hit != hit_end );

					if ( common_edge == ring.end() ) {
						// no sommon edge found, skip to the next
						continue;
					}
					// not part of the boundary
					ring.erase( common_edge );
					
					hit_end = hit;
					++hit;
				}
				
				// second pass: add each edge but the shared edge to the ring
				do {
					const CGAL::Point_3<Kernel>& s = hit->opposite()->vertex()->point();
					const CGAL::Point_3<Kernel>& t = hit->vertex()->point();
					
					CGAL::Segment_3<Kernel> candidate_segment( s, t );
					ring.push( candidate_segment );
					++hit;
				} while ( hit != hit_end );

				// mark the facet as processed
				ait->second.merged = true;
			}
		}

		// for each common plane
		for ( RingPartMap::const_iterator it = ring_parts.begin(); it != ring_parts.end(); ++it ) {
			LineString boundary;

			for ( Ring::const_iterator rit = it->second.begin(); rit != it->second.end(); ++rit ) {
				// filter out collinear points
				Ring::const_iterator next;
				next = rit;
				next++;
				if ( next == it->second.end() ) {
					next = it->second.begin();
				}

				if ( CGAL::collinear( rit->source(), rit->target(), next->target() ) ) {
						continue;
				}
				boundary.addPoint( rit->target() );
			}
			Polygon poly( boundary );
			addPolygon( poly );
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



