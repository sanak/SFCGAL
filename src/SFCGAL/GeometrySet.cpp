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

#include <SFCGAL/GeometrySet.h>

#include <SFCGAL/GeometryCollection.h>
#include <SFCGAL/Point.h>
#include <SFCGAL/MultiPoint.h>
#include <SFCGAL/LineString.h>
#include <SFCGAL/MultiLineString.h>
#include <SFCGAL/Triangle.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/MultiPolygon.h>
#include <SFCGAL/TriangulatedSurface.h>
#include <SFCGAL/PolyhedralSurface.h>
#include <SFCGAL/Solid.h>
#include <SFCGAL/MultiSolid.h>
#include <SFCGAL/Envelope.h>
#include <SFCGAL/TypeForDimension.h>

#include <SFCGAL/algorithm/covers.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

bool operator< ( const CGAL::Segment_2<SFCGAL::Kernel>& sega, const CGAL::Segment_2<SFCGAL::Kernel>& segb )
{
	if ( sega.source() == segb.source() ) {
		return sega.target() < segb.target();
	}
	return sega.source() < segb.source();
}

bool operator< ( const CGAL::Segment_3<SFCGAL::Kernel>& sega, const CGAL::Segment_3<SFCGAL::Kernel>& segb )
{
	if ( sega.source() == segb.source() ) {
		return sega.target() < segb.target();
	}
	return sega.source() < segb.source();
}

namespace SFCGAL {

	void _decompose_triangle( const Triangle& tri, typename GeometrySet<2>::SurfaceCollection& surfaces, dim_t<2> )
	{
		CGAL::Polygon_2<Kernel> outer;
		outer.push_back( tri.vertex(0).toPoint_2() );
		outer.push_back( tri.vertex(1).toPoint_2() );
		outer.push_back( tri.vertex(2).toPoint_2() );
		if ( outer.orientation() == CGAL::CLOCKWISE ) {
			outer.reverse_orientation();
		}
		surfaces.push_back( CGAL::Polygon_with_holes_2<Kernel>( outer ) );
	}
	void _decompose_triangle( const Triangle& tri, typename GeometrySet<3>::SurfaceCollection& surfaces, dim_t<3> )
	{
		CGAL::Triangle_3<Kernel> outtri( tri.vertex(0).toPoint_3(),
						 tri.vertex(1).toPoint_3(),
						 tri.vertex(2).toPoint_3() );
		surfaces.push_back( outtri );
	}

	void _decompose_polygon( const Polygon& poly, typename GeometrySet<2>::SurfaceCollection& surfaces, dim_t<2> )
	{
		surfaces.push_back( poly.toPolygon_with_holes_2() );
	}
	void _decompose_polygon( const Polygon& poly, typename GeometrySet<3>::SurfaceCollection& surfaces, dim_t<3> )
	{
		TriangulatedSurface surf;
		triangulate::triangulatePolygon3D( poly, surf );
		for ( size_t i = 0; i < surf.numTriangles(); ++i ) {
			const Triangle& tri = surf.triangleN( i );
			surfaces.push_back( CGAL::Triangle_3<Kernel>( tri.vertex(0).toPoint_3(),
								      tri.vertex(1).toPoint_3(),
									      tri.vertex(2).toPoint_3() )
					    );
		}
	}

	void _decompose_solid( const Solid& solid, typename GeometrySet<2>::VolumeCollection& volumes, dim_t<2> )
	{
		// nothing
	}
	void _decompose_solid( const Solid& solid, typename GeometrySet<3>::VolumeCollection& volumes, dim_t<3> )
	{
		volumes.push_back( *solid.exteriorShell().toPolyhedron_3<Kernel, MarkedPolyhedron >() );
	}

	template <int Dim>
	GeometrySet<Dim>::GeometrySet( )
	{
	}

	template <int Dim>
	GeometrySet<Dim>::GeometrySet( const Geometry& g )
	{
		_decompose( g );
	}

	template <int Dim>
	GeometrySet<Dim>::GeometrySet( const typename PrimitivePoint_d<Dim>::Type& g )
	{
		addPrimitive( g );
	}

	template <int Dim>
	GeometrySet<Dim>::GeometrySet( const typename PrimitiveSegment_d<Dim>::Type& g )
	{
		addPrimitive( g );
	}

	template <int Dim>
	GeometrySet<Dim>::GeometrySet( const typename PrimitiveSurface_d<Dim>::Type& g )
	{
		addPrimitive( g );
	}

	template <int Dim>
	GeometrySet<Dim>::GeometrySet( const typename PrimitiveVolume_d<Dim>::Type& g )
	{
		addPrimitive( g );
	}

	template <int Dim>
	void GeometrySet<Dim>::addGeometry( const Geometry& g )
	{
		_decompose( g );
	}

	template <int Dim>
	void GeometrySet<Dim>::addPrimitive( const PrimitiveBase<Dim>& p )
	{
		switch ( p.getType() )
		{
		case PrimitivePoint:
			_points.insert( p.template as<typename Point_d<Dim>::Type>() );
			break;
		case PrimitiveSegment:
			_segments.insert( p.template as<typename Segment_d<Dim>::Type>() );
			break;
		case PrimitiveSurface:
			_surfaces.push_back( p.template as<typename Surface_d<Dim>::Type>() );
			break;
		case PrimitiveVolume:
			_volumes.push_back( p.template as<typename Volume_d<Dim>::Type>() );
			break;
		}
	}

	template <>
	void GeometrySet<3>::addPrimitive( const CGAL::Object& o, bool pointsAsRing )
	{
		typedef typename Point_d<3>::Type TPoint;
		typedef typename Segment_d<3>::Type TSegment;
		typedef typename Surface_d<3>::Type TSurface;
		typedef typename Volume_d<3>::Type TVolume;
		if ( const TPoint * p = CGAL::object_cast<TPoint>( &o ) ) {
			_points.insert( TPoint( *p ) );
		}
		else if ( const std::vector<TPoint> * pts = CGAL::object_cast<std::vector<TPoint> >( &o ) ) {
			if ( pointsAsRing ) {
				// if pointsAsRing is true, build a polygon out of points
				// FIXME : we use triangulation here, which is not needed
				// We should have created a (planar) Polyhedron directly out of the points
				LineString ls;
				for ( size_t i = 0; i < pts->size(); ++i ) {
					ls.addPoint( (*pts)[i] );
				}
				Polygon poly(ls);
				_decompose_polygon( poly, _surfaces, dim_t<3>() );
			}
			else {
				std::copy( pts->begin(), pts->end(), std::inserter( _points, _points.end() ) );
			}
	        }
		else if ( const TSegment * p = CGAL::object_cast<TSegment>( &o ) ) {
			_segments.insert( TSegment( *p ) );
		}
		else if ( const TSurface * p = CGAL::object_cast<TSurface>( &o ) ) {
			_surfaces.push_back( TSurface( *p ) );
		}
		else if ( const TVolume * p = CGAL::object_cast<TVolume>( &o ) ) {
			_volumes.push_back( TVolume( *p ) );
		}
	}

	template <>
	void GeometrySet<2>::addPrimitive( const CGAL::Object& o, bool pointsAsRing )
	{
		typedef typename Point_d<2>::Type TPoint;
		typedef typename Segment_d<2>::Type TSegment;
		typedef typename Surface_d<2>::Type TSurface;
		typedef typename Volume_d<2>::Type TVolume;
		if ( const TPoint * p = CGAL::object_cast<TPoint>( &o ) ) {
			_points.insert( TPoint( *p ) );
		}
		else if ( const std::vector<TPoint> * pts = CGAL::object_cast<std::vector<TPoint> >( &o ) ) {
			if ( pointsAsRing ) {
				// if pointsAsRing is true, build a polygon out of points
				CGAL::Polygon_2<Kernel> poly;
				for ( size_t i = 0; i < pts->size(); ++i ) {
					poly.push_back( (*pts)[i] );
				}
				CGAL::Polygon_with_holes_2<Kernel> polyh( poly );
				_surfaces.push_back( polyh );
			}
			else {
				std::copy( pts->begin(), pts->end(), std::inserter( _points, _points.end() ) );
			}
	        }
		else if ( const CGAL::Triangle_2<Kernel>* tri = CGAL::object_cast<CGAL::Triangle_2<Kernel> >( &o ) ) {
			// convert to a polygon
			CGAL::Polygon_2<Kernel> poly;
			poly.push_back( tri->vertex(0) );
			poly.push_back( tri->vertex(1) );
			poly.push_back( tri->vertex(2) );
			CGAL::Polygon_with_holes_2<Kernel> polyh( poly );
			_surfaces.push_back( polyh );
		}
		else if ( const TSegment * p = CGAL::object_cast<TSegment>( &o ) ) {
			_segments.insert( TSegment( *p ) );
		}
		else if ( const TSurface * p = CGAL::object_cast<TSurface>( &o ) ) {
			_surfaces.push_back( TSurface( *p ) );
		}
		else if ( const TVolume * p = CGAL::object_cast<TVolume>( &o ) ) {
			_volumes.push_back( TVolume( *p ) );
		}
	}

	template <int Dim>
	void GeometrySet<Dim>::addPrimitive( const typename PrimitivePoint_d<Dim>::Type& p )
	{
		_points.insert( p );
	}

	template <int Dim>
	void GeometrySet<Dim>::addPrimitive( const typename PrimitiveSegment_d<Dim>::Type& p )
	{
		_segments.insert( p );
	}

	template <int Dim>
	void GeometrySet<Dim>::addPrimitive( const typename PrimitiveSurface_d<Dim>::Type& p )
	{
		_surfaces.push_back( p );
	}

	template <int Dim>
	void GeometrySet<Dim>::addPrimitive( const typename PrimitiveVolume_d<Dim>::Type& p )
	{
		_volumes.push_back( p );
	}

	template <int Dim>
	void GeometrySet<Dim>::_decompose( const Geometry& g )
	{
		if ( g.isEmpty() ) {
			return;
		}
		if ( g.is<GeometryCollection>() ) {
			const GeometryCollection& collect = g.as<GeometryCollection>();
			for ( size_t i = 0; i < g.numGeometries(); ++i ) {
				_decompose( collect.geometryN(i) );
			}
			return;
		}
		switch ( g.geometryTypeId() ) {
		case TYPE_POINT:
			_points.insert( g.as<Point>().toPoint_d<Dim>() );
			break;
		case TYPE_LINESTRING: {
			const LineString& ls = g.as<LineString>();
			for ( size_t i = 0; i < ls.numPoints() - 1; ++i ) {
				typename Segment_d<Dim>::Type seg( ls.pointN(i).toPoint_d<Dim>(),
									     ls.pointN(i+1).toPoint_d<Dim>() );
				_segments.insert( seg );
			}
			break;
		}
		case TYPE_TRIANGLE: {
			_decompose_triangle( g.as<Triangle>(), _surfaces, dim_t<Dim>() );
			break;
		}
		case TYPE_POLYGON: {
			_decompose_polygon( g.as<Polygon>(), _surfaces, dim_t<Dim>() );
			break;
		}
		case TYPE_TRIANGULATEDSURFACE: {
			const TriangulatedSurface& tri = g.as<TriangulatedSurface>();
			for ( size_t i = 0; i < tri.numTriangles(); ++i ) {
				_decompose( tri.triangleN( i ) );
			}
			break;
		}
		case TYPE_POLYHEDRALSURFACE: {
			const PolyhedralSurface& tri = g.as<PolyhedralSurface>();
			for ( size_t i = 0; i < tri.numPolygons(); ++i ) {
				_decompose( tri.polygonN( i ) );
			}
			break;
		}
		case TYPE_SOLID: {
			const Solid& solid = g.as<Solid>();
			_decompose_solid( solid, _volumes, dim_t<Dim>() );
			break;
		}
		default:
			break;
		}
	}

	// bbox of a 'volume' for 2D, will never be called
	CGAL::Bbox_2 computeBbox( const PrimitiveVolume_d<2>::Type& )
	{
		return CGAL::Bbox_2();
	}

	CGAL::Bbox_3 computeBbox( const PrimitiveVolume_d<3>::Type& p )
	{
		CGAL::Bbox_3 ret;
		MarkedPolyhedron::Point_const_iterator pit;
		const MarkedPolyhedron& vol = p.primitive();
		for ( pit = vol.points_begin(); pit != vol.points_end(); ++pit ) {
			ret = ret + pit->bbox();
		}
		return ret;
	}

	template <int Dim>
	void GeometrySet<Dim>::computeBoundingBoxes( typename HandleCollection<Dim>::Type& handles,
						     typename BoxCollection<Dim>::Type& boxes ) const
	{
		boxes.clear();
#if 0
		for ( typename PointCollection::const_iterator it = _points.begin(); it != _points.end(); ++it ) {
			handles.push_back( (PrimitiveBase<Dim>*)&*it );
			boxes.push_back( typename PrimitiveBox<Dim>::Type( it->primitive().bbox(), handles.back() ) );
		}
		for ( typename SegmentCollection::const_iterator it = _segments.begin(); it != _segments.end(); ++it ) {
			handles.push_back( (PrimitiveBase<Dim>*)&*it );
			boxes.push_back( typename PrimitiveBox<Dim>::Type( it->primitive().bbox(), handles.back() ) );
		}
		for ( typename SurfaceCollection::const_iterator it = _surfaces.begin(); it != _surfaces.end(); ++it ) {
			handles.push_back( (PrimitiveBase<Dim>*)&*it );
			boxes.push_back( typename PrimitiveBox<Dim>::Type( it->primitive().bbox(), handles.back() ) );
		}
		for ( typename VolumeCollection::const_iterator it = _volumes.begin(); it != _volumes.end(); ++it ) {
			handles.push_back( (PrimitiveBase<Dim>*)&*it );
			boxes.push_back( typename PrimitiveBox<Dim>::Type( computeBbox( *it ),
									   handles.back() ) );
		}
#endif
		for ( ConstIterator it = primitives_begin(); it != primitives_end(); ++it ) {
			handles.push_back( &*it );
			boxes.push_back( typename PrimitiveBox<Dim>::Type( it->bbox(), handles.back() ) );
		}
	}

	template <int Dim>
	void recompose_points( const typename GeometrySet<Dim>::PointCollection& points,
			       std::vector<Geometry *>& rpoints,
			       dim_t<Dim> )
	{
		if ( points.empty() ) {
			return;
		}
		else {
			for ( typename GeometrySet<Dim>::PointCollection::const_iterator it = points.begin();
			      it != points.end();
			      ++it ) {
				rpoints.push_back( new Point( it->primitive() ) );
			}
		}
	}


	template <int Dim>
	void recompose_segments( const typename GeometrySet<Dim>::SegmentCollection& segments,
				 std::vector<Geometry*>& lines,
				 dim_t<Dim> )
	{
		if ( segments.empty() ) {
			return;
		}

		// convert segments to linestring / multi linestring
		LineString* ls = new LineString;
		bool first = true;
		typename Segment_d<Dim>::Type lastSeg;
		for ( typename GeometrySet<Dim>::SegmentCollection::const_iterator it = segments.begin();
		      it != segments.end();
		      ++it ) {
			if ( !first && lastSeg.target() != it->primitive().source() ) {
				lines.push_back( ls );
				ls = new LineString;
				first = true;
			}
			if ( first ) {
				ls->addPoint( new Point( it->primitive().source() ) );
				first = false;
			}
			ls->addPoint( new Point( it->primitive().target() ) );
			lastSeg = it->primitive();
		}
		lines.push_back( ls );
	}

	void recompose_surfaces( const GeometrySet<2>::SurfaceCollection& surfaces, std::vector<Geometry*>& output, dim_t<2> )
	{
		for ( GeometrySet<2>::SurfaceCollection::const_iterator it = surfaces.begin(); it != surfaces.end(); ++it ) {
			if ( it->primitive().holes_begin() == it->primitive().holes_end() &&
			     it->primitive().outer_boundary().size() == 3 ) {
				CGAL::Polygon_2<Kernel>::Vertex_iterator vit = it->primitive().outer_boundary().vertices_begin();
				CGAL::Point_2<Kernel> p1( *vit++ );
				CGAL::Point_2<Kernel> p2( *vit++ );
				CGAL::Point_2<Kernel> p3( *vit++ );
				output.push_back( new Triangle(CGAL::Triangle_2<Kernel>( p1, p2, p3 )) );
			}
			else {
				output.push_back( new Polygon( it->primitive() ) );
			}
		}
	}

	void recompose_surfaces( const GeometrySet<3>::SurfaceCollection& surfaces, std::vector<Geometry*>& output, dim_t<3> )
	{
		if ( surfaces.empty() ) {
			return;
		}
		// TODO : regroup triangles of the same mesh
		if ( surfaces.size() == 1 ) {
			output.push_back( new Triangle( surfaces.begin()->primitive() ) );
			return;
		}
		TriangulatedSurface *tri = new TriangulatedSurface;
		for ( GeometrySet<3>::SurfaceCollection::const_iterator it = surfaces.begin(); it != surfaces.end(); ++it ) {
			tri->addTriangle( new Triangle( it->primitive() ) );
		}
		output.push_back( tri );
	}

	void recompose_volumes( const GeometrySet<2>::VolumeCollection& volumes, std::vector<Geometry*>& output, dim_t<2> )
	{
	}

	void recompose_volumes( const GeometrySet<3>::VolumeCollection& volumes, std::vector<Geometry*>& output, dim_t<3> )
	{
		if ( volumes.empty() ) {
			return;
		}
		for ( GeometrySet<3>::VolumeCollection::const_iterator vit = volumes.begin(); vit != volumes.end(); ++vit ) {
			if ( vit->flags() & FLAG_IS_PLANAR ) {
				// extract the boundary
				std::list<CGAL::Point_3<Kernel> > boundary;
				for ( MarkedPolyhedron::Halfedge_const_iterator it = vit->primitive().halfedges_begin();
				      it != vit->primitive().halfedges_end();
				      ++it ) {
					if ( !it->is_border() )
						continue;
					CGAL::Point_3<Kernel> p1 = it->prev()->vertex()->point();
					CGAL::Point_3<Kernel> p2 = it->vertex()->point();
					// TODO: test for colinearity
					// Temporary vertice may have been introduced during triangulations
					// and since we expect here a planar surface, it is safe to simplify the boundary
					// by eliminating collinear points.
					if ( boundary.size() == 0 ) {
						boundary.push_back( p1 );
						boundary.push_back( p2 );
					}
					else if ( boundary.back() == p1 ) {
						boundary.push_back( p2 );
					}
					else if ( boundary.front() == p2 ) {
						boundary.push_front( p1 );
					}
				}
				if ( boundary.size() == 3 ) {
					// It is a triangle
					
					Point p[3];
					std::list<CGAL::Point_3<Kernel> >::const_iterator it = boundary.begin();
					for ( size_t i = 0; i < 3; ++i, ++it ) {
						p[i] = *it;
					}
					output.push_back( new Triangle( p[0], p[1], p[2] ) );
				}
				else {
					// Else it is a polygon
					LineString *ls = new LineString;
					for ( std::list<CGAL::Point_3<Kernel> >::const_iterator it = boundary.begin(); it != boundary.end(); ++it ) {
						ls->addPoint( *it );
					}
					output.push_back( new Polygon( ls ) );
				}
			}
			else {
				
				PolyhedralSurface *shell = new PolyhedralSurface( vit->primitive() );
				// TODO: test open / closed
				output.push_back( new Solid( shell ) );
			}
		}
	}

	template <int Dim>
	std::auto_ptr<Geometry> GeometrySet<Dim>::recompose() const
	{
		std::vector<Geometry*> geometries;

		recompose_points( _points, geometries, dim_t<Dim>() );
		recompose_segments( _segments, geometries, dim_t<Dim>() );
		recompose_surfaces( _surfaces, geometries, dim_t<Dim>() );
		recompose_volumes( _volumes, geometries, dim_t<Dim>() );

		if ( geometries.empty() ) {
			return std::auto_ptr<Geometry>( new GeometryCollection );
		}
		if ( geometries.size() == 1 ) {
			return std::auto_ptr<Geometry>( geometries[0] );
		}
		// else we have a mix of different types
		bool hasCommonType = true;
		int commonType = geometries[0]->geometryTypeId();
		for ( size_t i = 0; i < geometries.size(); ++i ) {
			if ( geometries[i]->geometryTypeId() != commonType ) {
				hasCommonType = false;
				break;
			}
		}
		GeometryCollection* ret;
		if ( hasCommonType ) {
			if ( commonType == TYPE_POINT ) {
				ret = new MultiPoint;
			}
			else if ( commonType == TYPE_LINESTRING ) {
				ret = new MultiLineString;
			}
			else if ( commonType == TYPE_POLYGON ) {
				ret = new MultiPolygon;
			}
			else if ( commonType == TYPE_SOLID ) {
				ret = new MultiSolid;
			}
		}
		else {
			ret = new GeometryCollection;
		}
		for ( size_t i = 0; i < geometries.size(); ++i ) {
			ret->addGeometry( geometries[i] );
		}
		return std::auto_ptr<Geometry>( ret );
	}

	void _collect_points( const CGAL::Polygon_with_holes_2<Kernel>& poly, GeometrySet<2>::PointCollection& points )
	{
		for ( CGAL::Polygon_2<Kernel>::Vertex_iterator vit = poly.outer_boundary().vertices_begin();
		      vit != poly.outer_boundary().vertices_end();
		      ++vit ) {
			points.insert( *vit );
		}
		for ( CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator hit = poly.holes_begin();
		      hit != poly.holes_end();
		      ++hit ) {
			for ( CGAL::Polygon_2<Kernel>::Vertex_iterator vit = hit->vertices_begin();
			      vit != hit->vertices_end();
			      ++vit ) {
				points.insert( *vit );
			}
		}
	}

	void _collect_points( const CGAL::Triangle_3<Kernel>& tri, GeometrySet<3>::PointCollection& points )
	{
		points.insert( tri.vertex(0) );
		points.insert( tri.vertex(1) );
		points.insert( tri.vertex(2) );
	}

	void _collect_points( const NoVolume&, GeometrySet<2>::PointCollection& points )
	{
		// nothing
	}

	void _collect_points( const MarkedPolyhedron& poly, GeometrySet<3>::PointCollection& points )
	{
		for ( MarkedPolyhedron::Vertex_const_iterator vit = poly.vertices_begin();
		      vit != poly.vertices_end();
		      ++vit ) {
			points.insert( vit->point() );
		}
	}

	template <int Dim>
	void GeometrySet<Dim>::collectPoints( const PrimitiveBase<Dim>& pa )
	{
		typedef typename Point_d<Dim>::Type TPoint;
		typedef typename Segment_d<Dim>::Type TSegment;
		typedef typename Surface_d<Dim>::Type TSurface;
		typedef typename Volume_d<Dim>::Type TVolume;

		switch ( pa.getType() ) {
		case PrimitivePoint: {
			_points.insert( pa.template as<TPoint>() );
			break;
		}
		case PrimitiveSegment: {
			const TSegment& seg = pa.template as<TSegment>();
			_points.insert( seg.source() );
			_points.insert( seg.target() );
			break;
		}
		case PrimitiveSurface: {
			_collect_points( pa.template as<TSurface>(), _points );
			break;
		}
		case PrimitiveVolume: {
			_collect_points( pa.template as<TVolume>(), _points );
			break;
		}
		}
	}

	template <int Dim, class IT>
	void _filter_covered( IT ibegin, IT iend, GeometrySet<Dim>& output )
	{
		for ( IT it = ibegin; it != iend; ++it ) {
			GeometrySet<Dim> v1;
			v1.addPrimitive( *it );
			bool v1_covered = false;
			for ( IT it2 = it; it2 != iend; ++it2 ) {
				if ( it == it2 ) {
					continue;
				}

				GeometrySet<Dim> v2;
				v2.addPrimitive( *it2 );
				if ( algorithm::covers( v2, v1 ) ) {
					v1_covered = true;
					break;
				}
			}
			// if its not covered by another volume
			if ( !v1_covered ) {
				// and not covered by another already inserted primitive
				bool b = algorithm::covers( output, v1 );
				if ( !b ) {
					output.addPrimitive( *it );
				}
			}
		}
	}

	template <int Dim>
	void GeometrySet<Dim>::filterCovered( GeometrySet<Dim>& output ) const
	{
		_filter_covered( _volumes.begin(), _volumes.end(), output );
		_filter_covered( _surfaces.begin(), _surfaces.end(), output );
		_filter_covered( _segments.begin(), _segments.end(), output );
		_filter_covered( _points.begin(), _points.end(), output );
	}

	std::ostream& operator<<( std::ostream& ostr, const GeometrySet<2>& g )
	{
		ostr << "points: ";
 		std::ostream_iterator<PrimitivePoint_d<2>::Type > out_pt (ostr,", ");
		std::copy( g.points().begin(), g.points().end(), out_pt );
		ostr << std::endl << "segments: ";
		std::ostream_iterator<PrimitiveSegment_d<2>::Type > out_seg (ostr,", ");
		std::copy( g.segments().begin(), g.segments().end(), out_seg );
		ostr << std::endl << "surfaces: ";
		std::ostream_iterator<PrimitiveSurface_d<2>::Type > out_surf (ostr,", ");
		std::copy( g.surfaces().begin(), g.surfaces().end(), out_surf );
		ostr << std::endl;
		return ostr;
	}

	std::ostream& operator<<( std::ostream& ostr, const GeometrySet<3>& g )
	{
		ostr << "points: ";
 		std::ostream_iterator<PrimitivePoint_d<3>::Type > out_pt (ostr,", ");
		std::copy( g.points().begin(), g.points().end(), out_pt );
		ostr << std::endl << "segments: ";
		std::ostream_iterator<PrimitiveSegment_d<3>::Type > out_seg (ostr,", ");
		std::copy( g.segments().begin(), g.segments().end(), out_seg );
		ostr << std::endl << "surfaces: ";
		std::ostream_iterator<PrimitiveSurface_d<3>::Type > out_surf (ostr,", ");
		std::copy( g.surfaces().begin(), g.surfaces().end(), out_surf );
		ostr << std::endl << "volumes: ";
		std::ostream_iterator<PrimitiveVolume_d<3>::Type > out_vol (ostr,", ");
		std::copy( g.volumes().begin(), g.volumes().end(), out_vol );
		ostr << std::endl;
		return ostr;
	}

	template <int Dim>
	template <class T>
	void GeometrySet<Dim>::IteratorBase<T>::increment()
	{
		if ( pointIt != pointEnd ) {
			++pointIt;
			return;
		}
		if ( segmentIt != segmentEnd ) {
			++segmentIt;
			return;
		}
		if ( surfaceIt != surfaceEnd ) {
			++surfaceIt;
			return;
		}
		if ( volumeIt != volumeEnd ) {
			++volumeIt;
		}
	}

	template <int Dim>
	template <class T>
	GeometrySet<Dim>::IteratorBase<T>::IteratorBase( PointIterator pit,
							 PointIterator pend,
							 SegmentIterator lit,
							 SegmentIterator lend,
							 SurfaceIterator sit,
							 SurfaceIterator send,
							 VolumeIterator vit,
							 VolumeIterator vend ) :
		pointIt( pit ), pointEnd( pend ),
		segmentIt( lit ), segmentEnd( lend ),
		surfaceIt( sit ), surfaceEnd( send ),
		volumeIt( vit ), volumeEnd( vend )
	{}

	template <int Dim>
	template <class T>
	bool GeometrySet<Dim>::IteratorBase<T>::equal( const IteratorBase<T>& other ) const
	{
		return pointIt == other.pointIt &&
			segmentIt == other.segmentIt &&
			surfaceIt == other.surfaceIt &&
			volumeIt == other.volumeIt &&
			pointEnd == other.pointEnd &&
			segmentEnd == other.segmentEnd &&
			surfaceEnd == other.surfaceEnd &&
			volumeEnd == other.volumeEnd;
	}
	
	template <int Dim>
	template <class T>
	T& GeometrySet<Dim>::IteratorBase<T>::dereference() const
	{
		if ( pointIt != pointEnd ) {
			return (T&)*pointIt;
		}
		if ( segmentIt != segmentEnd ) {
			return (T&)*segmentIt;
		}
		if ( surfaceIt != surfaceEnd ) {
			return (T&)*surfaceIt;
		}
		// if ( volumeIt != volumeEnd ) {
		return (T&)*volumeIt;
	}

	template <int Dim>
	typename GeometrySet<Dim>::ConstIterator GeometrySet<Dim>::primitives_begin( int selection ) const
	{
		if ( selection == -1 ) {
			return ConstIterator( _points.begin(), _points.end(),
					      _segments.begin(), _segments.end(),
					      _surfaces.begin(), _surfaces.end(),
					      _volumes.begin(), _volumes.end() );
		}
		if ( selection == PrimitivePoint ) {
			return ConstIterator( _points.begin(), _points.end(),
					      _segments.end(), _segments.end(),
					      _surfaces.end(), _surfaces.end(),
					      _volumes.end(), _volumes.end() );
		}
		if ( selection == PrimitiveSegment ) {
			return ConstIterator( _points.end(), _points.end(),
					      _segments.begin(), _segments.end(),
					      _surfaces.end(), _surfaces.end(),
					      _volumes.end(), _volumes.end() );
		}
		if ( selection == PrimitiveSurface ) {
			return ConstIterator( _points.end(), _points.end(),
					      _segments.end(), _segments.end(),
					      _surfaces.begin(), _surfaces.end(),
					      _volumes.end(), _volumes.end() );
		}
		// if ( selection == PrimitiveVolume ) {
		return ConstIterator( _points.end(), _points.end(),
				      _segments.end(), _segments.end(),
				      _surfaces.end(), _surfaces.end(),
				      _volumes.begin(), _volumes.end() );
	}

	template <int Dim>
	typename GeometrySet<Dim>::Iterator GeometrySet<Dim>::primitives_begin( int selection )
	{
		if ( selection == -1 ) {
			return Iterator( _points.begin(), _points.end(),
					      _segments.begin(), _segments.end(),
					      _surfaces.begin(), _surfaces.end(),
					      _volumes.begin(), _volumes.end() );
		}
		if ( selection == PrimitivePoint ) {
			return Iterator( _points.begin(), _points.end(),
					      _segments.end(), _segments.end(),
					      _surfaces.end(), _surfaces.end(),
					      _volumes.end(), _volumes.end() );
		}
		if ( selection == PrimitiveSegment ) {
			return Iterator( _points.end(), _points.end(),
					      _segments.begin(), _segments.end(),
					      _surfaces.end(), _surfaces.end(),
					      _volumes.end(), _volumes.end() );
		}
		if ( selection == PrimitiveSurface ) {
			return Iterator( _points.end(), _points.end(),
					      _segments.end(), _segments.end(),
					      _surfaces.begin(), _surfaces.end(),
					      _volumes.end(), _volumes.end() );
		}
		//		if ( selection == PrimitiveVolume ) {
		return Iterator( _points.end(), _points.end(),
				 _segments.end(), _segments.end(),
				 _surfaces.end(), _surfaces.end(),
				 _volumes.begin(), _volumes.end() );
	}

	template <int Dim>
	typename GeometrySet<Dim>::ConstIterator GeometrySet<Dim>::primitives_end() const
	{
		return ConstIterator( _points.end(), _points.end(),
				      _segments.end(), _segments.end(),
				      _surfaces.end(), _surfaces.end(),
				      _volumes.end(), _volumes.end() );
	}
	
	template <int Dim>
	typename GeometrySet<Dim>::Iterator GeometrySet<Dim>::primitives_end()
	{
		return Iterator( _points.end(), _points.end(),
				 _segments.end(), _segments.end(),
				 _surfaces.end(), _surfaces.end(),
				 _volumes.end(), _volumes.end() );
	}
	
	template <int Dim>
	typename GeometrySet<Dim>::Iterator GeometrySet<Dim>::insert( typename GeometrySet<Dim>::Iterator it, const PrimitiveBase<Dim>& value )
	{
		switch ( value.getType() ) {
		case PrimitivePoint:
			_points.insert( it.pointIt, value.template as<typename Point_d<Dim>::Type>() );
			break;
		case PrimitiveSegment:
			_segments.insert( it.segmentIt, value.template as<typename Segment_d<Dim>::Type>() );
			break;
		case PrimitiveSurface:
			_surfaces.insert( it.surfaceIt, value.template as<typename Surface_d<Dim>::Type>() );
			break;
		case PrimitiveVolume:
			_volumes.insert( it.volumeIt, value.template as<typename Volume_d<Dim>::Type>() );
			break;
		}
		return it;
	}

	template <int Dim>
	typename GeometrySet<Dim>::ConstIterator GeometrySet<Dim>::insert( typename GeometrySet<Dim>::ConstIterator it, const PrimitiveBase<Dim>& value )
	{
		return it;
	}

	template class GeometrySet<2>;
	template class GeometrySet<3>;
} // SFCGAL
