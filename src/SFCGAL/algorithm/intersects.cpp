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
#include <map>
#include <sstream>

#include <SFCGAL/Kernel.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/algorithm/covers.h>
#include <SFCGAL/triangulate/triangulateInGeometrySet.h>
#include <SFCGAL/GeometrySet.h>
#include <SFCGAL/Envelope.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/box_intersection_d.h>

//#define CACHE_TRIANGULATION

namespace SFCGAL {
namespace algorithm
{
	//
	// Type of pa must be of larger dimension than type of pb
	bool _intersects( const PrimitiveBase& pa, const PrimitiveBase& pb, const dim_t<2>& )
	{
		//
		// Point vs. Point
		//

		if ( pa.getType() == PrimitivePoint && pb.getType() == PrimitivePoint ) {
			return pa.as<PrimitivePoint_d<2>::Type>().primitive()
				== pb.as<PrimitivePoint_d<2>::Type>().primitive();
		}

		//
		// Segment vs. Point
		//

		else if ( pa.getType() == PrimitiveSegment && pb.getType() == PrimitivePoint ) {
			const CGAL::Segment_2<Kernel>& seg = pa.as<PrimitiveSegment_d<2>::Type>().primitive();
			const CGAL::Point_2<Kernel>& pt = pb.as<PrimitivePoint_d<2>::Type>().primitive();
			return seg.has_on( pt );
		}

		//
		// Segment vs. Segment
		//

		else if ( pa.getType() == PrimitiveSegment && pb.getType() == PrimitiveSegment ) {
			const CGAL::Segment_2<Kernel>& seg1 = pa.as<PrimitiveSegment_d<2>::Type>().primitive();
			const CGAL::Segment_2<Kernel>& seg2 = pb.as<PrimitiveSegment_d<2>::Type>().primitive();
			return CGAL::do_intersect( seg1, seg2 );
		}
		//
		// Polygon vs. Point
		//

		else if ( pa.getType() == PrimitiveSurface && pb.getType() == PrimitivePoint ) {
			// Polygon versus Point
			const CGAL::Polygon_with_holes_2<Kernel> *poly = &pa.as<PrimitiveSurface_d<2>::Type>().primitive();
			const CGAL::Point_2<Kernel> *pt = &pb.as<PrimitivePoint_d<2>::Type>().primitive();

			int b1 = poly->outer_boundary().bounded_side( *pt );
			if ( b1 == CGAL::ON_BOUNDARY ) {
				return true;
			}
			if ( b1 == CGAL::ON_BOUNDED_SIDE ) {
				CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it;
				for ( it = poly->holes_begin(); it != poly->holes_end(); ++it ) {
					int b = it->bounded_side( *pt );
					if ( b == CGAL::ON_BOUNDED_SIDE ) {
						return false;
					}
				}
			}
			else {
				return false;
			}
			return true;
		}

		//
		// Polygon vs. Segment
		//

		else if ( pa.getType() == PrimitiveSurface && pb.getType() == PrimitiveSegment ) {
			const CGAL::Polygon_with_holes_2<Kernel> *poly = &pa.as<PrimitiveSurface_d<2>::Type>().primitive();
			const CGAL::Segment_2<Kernel> *seg = &pb.as<PrimitiveSegment_d<2>::Type>().primitive();
			
			// 1. if the segment intersects a boundary of the polygon, returns true
			// 2. else, if one of the point of the segment intersects the polygon, returns true

			GeometrySet<2> segment;
			segment.addSegments( seg, seg+1 );

			// 1.
			GeometrySet<2> boundary;
			boundary.addSegments( poly->outer_boundary().edges_begin(),
			 		      poly->outer_boundary().edges_end() );

			// recurse call
			if ( intersects( boundary, segment ) ) {
			 	return true;
			}
			for ( CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it = poly->holes_begin();
			      it != poly->holes_end();
			      ++it ) {
				GeometrySet<2> hole;
				hole.addSegments( it->edges_begin(), it->edges_end() );

				// recurse call
				if ( intersects( hole, segment ) ) {
					return true;
				}
			}

			// 2. call the polygon, point version
			CGAL::Point_2<Kernel> pt = seg->source();
			PrimitivePoint_d<2>::Type ppt( pt );
			return intersects( pa, ppt );
		}

		//
		// Polygon vs. Polygon
		//

		else if ( pa.getType() == PrimitiveSurface && pb.getType() == PrimitiveSurface ) {
			const CGAL::Polygon_with_holes_2<Kernel> *poly1 = &pa.as<PrimitiveSurface_d<2>::Type>().primitive();
			const CGAL::Polygon_with_holes_2<Kernel> *poly2 = &pb.as<PrimitiveSurface_d<2>::Type>().primitive();

			// 1. if rings intersects, returns true
			// 2. else, if poly1 is inside poly2 or poly1 inside poly2 (but not in holes), returns true

			GeometrySet<2> rings1, rings2;
			rings1.addSegments( poly1->outer_boundary().edges_begin(),
					    poly1->outer_boundary().edges_end() );
			for ( CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it = poly1->holes_begin();
			      it != poly1->holes_end();
			      ++it ) {
				rings1.addSegments( it->edges_begin(), it->edges_end() );
			}
			rings2.addSegments( poly2->outer_boundary().edges_begin(),
					    poly2->outer_boundary().edges_end() );
			for ( CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it = poly2->holes_begin();
			      it != poly2->holes_end();
			      ++it ) {
				rings2.addSegments( it->edges_begin(), it->edges_end() );
			}

			// 1.
			if ( intersects( rings1, rings2 ) ) {
				return true;
			}

			// 2.
			CGAL::Bbox_2 box1, box2;
			box1 = poly1->bbox();
			box2 = poly2->bbox();
			Envelope e1( box1.xmin(), box1.xmax(), box1.ymin(), box1.ymax() );
			Envelope e2( box2.xmin(), box2.xmax(), box2.ymin(), box2.ymax() );

			// if pa is inside pb
			if ( Envelope::contains( e2, e1 ) ) {
				// is pa inside one of pb's holes ?
				CGAL::Point_2<Kernel> pt = *poly1->outer_boundary().vertices_begin();
				for ( CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it = poly2->holes_begin();
				      it != poly2->holes_end();
				      ++it ) {
					CGAL::Bounded_side b2 = it->bounded_side( pt );
					if ( b2 == CGAL::ON_BOUNDED_SIDE ) {
						return false;
					}
				}
				return true;
			}
			// if pb is inside pa
			if ( Envelope::contains( e1, e2 ) ) {
				// is pa inside one of pb's holes ?
				CGAL::Point_2<Kernel> pt = *poly2->outer_boundary().vertices_begin();
				for ( CGAL::Polygon_with_holes_2<Kernel>::Hole_const_iterator it = poly1->holes_begin();
				      it != poly1->holes_end();
				      ++it ) {
					CGAL::Bounded_side b2 = it->bounded_side( pt );
					if ( b2 == CGAL::ON_BOUNDED_SIDE ) {
						return false;
					}
				}
				return true;
			}
			return false;
		}
		return false;
	}


	//
	// intersects of a volume with any other type
	bool intersects_volume_x( const MarkedPolyhedron* polyhedron, const PrimitiveBase& geometry )
	{
		typedef CGAL::Polyhedral_mesh_domain_3<MarkedPolyhedron, Kernel> Mesh_domain;
		
		// intersection between a solid and a geometry
		// 1. either one of the geometry' point lies inside the solid
		// 2. or the geometry intersects one of the surfaces
		
		// 1.
		
		Mesh_domain ext_domain( *polyhedron );
		Mesh_domain::Is_in_domain is_in_poly( ext_domain );
		
		GeometrySet<3> points;
		points.collectPoints( geometry );
		for ( GeometrySet<3>::PointCollection::const_iterator pit = points.points().begin();
		      pit != points.points().end(); ++pit ) {
			if ( is_in_poly( pit->primitive() ) ) {
				return true;
			}
		}
		
		// 2.
		
		// triangulate the polyhedron_3
		GeometrySet<3> g;
		g.addPrimitive( geometry );
		
		GeometrySet<3> triangles;
		triangulate::triangulate( *polyhedron, triangles );
		
		return intersects( g, triangles );
	}

	//
	// Type of pa must be of larger dimension than type of pb
	bool _intersects( const PrimitiveBase& pa, const PrimitiveBase& pb, const dim_t<3>& )
	{
		if ( pa.getType() == PrimitivePoint && pb.getType() == PrimitivePoint ) {
			return pa.as<PrimitivePoint_d<3>::Type>().primitive()
				== pb.as<PrimitivePoint_d<3>::Type>().primitive();
		}
		else if ( pa.getType() == PrimitiveSegment && pb.getType() == PrimitivePoint ) {
			const CGAL::Segment_3<Kernel>& seg = pa.as<PrimitiveSegment_d<3>::Type>().primitive();
			const CGAL::Point_3<Kernel>& pt = pb.as<PrimitivePoint_d<3>::Type>().primitive();
			return seg.has_on( pt );
		}
		else if ( pa.getType() == PrimitiveSegment && pb.getType() == PrimitiveSegment ) {
			const CGAL::Segment_3<Kernel>& sega = pa.as<PrimitiveSegment_d<3>::Type>().primitive();
			const CGAL::Segment_3<Kernel>& segb = pb.as<PrimitiveSegment_d<3>::Type>().primitive();
			return CGAL::do_intersect( sega, segb );
		}
		if ( pa.getType() == PrimitiveVolume ) {
			return intersects_volume_x( &pa.as<PrimitiveVolume_d<3>::Type>().primitive(), pb );
		}
		if ( pa.getType() == PrimitiveSurface && pb.getType() == PrimitivePoint ) {
			const CGAL::Triangle_3<Kernel>& tri = pa.as<PrimitiveSurface_d<3>::Type>().primitive();
			const CGAL::Point_3<Kernel>& pt = pb.as<PrimitivePoint_d<3>::Type>().primitive();
			return tri.has_on( pt );
		}
		if ( pa.getType() == PrimitiveSurface && pb.getType() == PrimitiveSegment ) {
			const CGAL::Triangle_3<Kernel>& tri = pa.as<PrimitiveSurface_d<3>::Type>().primitive();
			const CGAL::Segment_3<Kernel>& seg = pb.as<PrimitiveSegment_d<3>::Type>().primitive();
			return CGAL::do_intersect( tri, seg );
		}
		if ( pa.getType() == PrimitiveSurface && pb.getType() == PrimitiveSurface ) {
			const CGAL::Triangle_3<Kernel>& tri1 = pa.as<PrimitiveSurface_d<3>::Type>().primitive();
			const CGAL::Triangle_3<Kernel>& tri2 = pb.as<PrimitiveSurface_d<3>::Type>().primitive();
			return CGAL::do_intersect( tri1, tri2 );
		}
		return false;
	}

	//
	// We deal here with symmetric call
	template <int Dim>
	bool dispatch_intersects_sym( const PrimitiveBase& pa, const PrimitiveBase& pb, const dim_t<Dim>& dim )
	{
		// assume types are ordered by dimension within the boost::variant
		if ( pa.getType() >= pb.getType() ) {
			return _intersects( pa, pb, dim );
		}
		else {
			return _intersects( pb, pa, dim );
		}
	}

	bool intersects( const PrimitiveBase& pa, const PrimitiveBase& pb )
	{
		if ( pa.is3D() ) {
			return dispatch_intersects_sym( pa, pb, dim_t<3>() );
		}
		return dispatch_intersects_sym( pa, pb, dim_t<2>() );
	}

	struct found_an_intersection{};

	template <int Dim>
	struct intersects_cb
	{
		void operator()( const typename PrimitiveBox<Dim>::Type& a,
				 const typename PrimitiveBox<Dim>::Type& b )
		{
			if ( dispatch_intersects_sym( *a.handle(), *b.handle(), dim_t<Dim>() ) ) {
				throw found_an_intersection();
			}
		}
	};

	template <int Dim>
	bool intersects( const GeometrySet<Dim>& a, const GeometrySet<Dim>& b )
	{
		typename SFCGAL::HandleCollection ahandles, bhandles;
		typename SFCGAL::BoxCollection<Dim>::Type aboxes, bboxes;
		a.computeBoundingBoxes( ahandles, aboxes );
		b.computeBoundingBoxes( bhandles, bboxes );

		try {
			intersects_cb<Dim> cb;
			CGAL::box_intersection_d( aboxes.begin(), aboxes.end(),
						  bboxes.begin(), bboxes.end(),
						  cb );
		}
		catch ( found_an_intersection& e ) {
			return true;
		}
		return false;
	}

	template bool intersects<2>( const GeometrySet<2>& a, const GeometrySet<2>& b );
	template bool intersects<3>( const GeometrySet<3>& a, const GeometrySet<3>& b );

	bool intersects( const Geometry& ga, const Geometry& gb )
	{
		GeometrySet<2> gsa( ga );
		GeometrySet<2> gsb( gb );

		return intersects( gsa, gsb );
	}

	bool intersects3D( const Geometry& ga, const Geometry& gb )
	{
		GeometrySet<3> gsa( ga );
		GeometrySet<3> gsb( gb );

		return intersects( gsa, gsb );
	}

}
}
