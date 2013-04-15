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

#include <CGAL/box_intersection_d.h>
#include <CGAL/corefinement_operations.h>

#include <SFCGAL/algorithm/difference.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/Geometry.h>
#include <SFCGAL/Polygon.h>
#include <SFCGAL/GeometrySet.h>

// #include <CGAL/Boolean_set_operations_2.h>

namespace SFCGAL {

//
// FIXME: workaround CGAL Boolean_set_operations_2 / intersections inclusion bug
//
typedef std::list<CGAL::Polygon_with_holes_2<Kernel> > PolygonList;
void CGAL_poly_difference( const CGAL::Polygon_2<Kernel> p1, const CGAL::Polygon_2<Kernel> p2, PolygonList& );
void CGAL_poly_difference( const CGAL::Polygon_with_holes_2<Kernel> p1, const CGAL::Polygon_with_holes_2<Kernel> p2, GeometrySet<2>::SurfaceCollection& );

namespace algorithm
{
	// get the point position on segment
	// in [0,1] if its on the interior
	static Kernel::FT point_position( const typename Segment_d<2>::Type& seg, const typename Point_d<2>::Type& pt )
	{
		const typename Point_d<2>::Type pA = seg.source();
		const typename Point_d<2>::Type pB = seg.target();
		Kernel::FT num = pt.x() - pA.x();
		Kernel::FT den = pB.x() - pA.x();
		if ( den == 0 ) {
			num = pt.y() - pA.y();
			den = pB.y() - pA.y();
		}
		return num / den;
	}
	static Kernel::FT point_position( const typename Segment_d<3>::Type& seg, const typename Point_d<3>::Type& pt )
	{
		const typename Point_d<3>::Type pA = seg.source();
		const typename Point_d<3>::Type pB = seg.target();
		Kernel::FT num = pt.x() - pA.x();
		Kernel::FT den = pB.x() - pA.x();
		if ( den == 0 ) {
			num = pt.y() - pA.y();
			den = pB.y() - pA.y();
		}
		return num / den;
	}

	// prim intersects seg
	template <int Dim>
	static void substract_from_segment( const PrimitiveBase<Dim>& seg, const PrimitiveBase<Dim>& prim, GeometrySet<Dim>& output )
	{
		// seg is a segment
		if ( prim.getType() == PrimitiveSegment ) {
			const typename Segment_d<Dim>::Type& sega = seg.template as<typename Segment_d<Dim>::Type>();
			const typename Segment_d<Dim>::Type& segb = prim.template as<typename Segment_d<Dim>::Type>();
			typename Point_d<Dim>::Type pA( sega.source() );
			typename Point_d<Dim>::Type pB( sega.target() );
			typename Point_d<Dim>::Type pC( segb.source() );
			typename Point_d<Dim>::Type pD( segb.target() );
			Kernel::FT sC = point_position( sega, pC );
			Kernel::FT sD = point_position( sega, pD );
			if ( sC > sD ) {
				std::swap( sC, sD );
				std::swap( pC, pD );
			}
			
			if ( sC > 0 ) {
				output.addPrimitive( typename Segment_d<Dim>::Type( pA, pC ) );
			}
			if ( sD < 1 ) {
				output.addPrimitive( typename Segment_d<Dim>::Type( pD, pB ) );
			}
		}
		else {
			output.addPrimitive( seg );
		}
	}

	// prim intersects surf
	static void substract_from_surface( const PrimitiveBase<2>& surf, const PrimitiveBase<2>& prim, GeometrySet<2>& output )
	{
		if ( prim.getType() == PrimitiveSurface ) {
			const Surface_d<2>::Type& poly1 = surf.as<Surface_d<2>::Type>();
			const Surface_d<2>::Type& poly2 = prim.as<Surface_d<2>::Type>();

			CGAL_poly_difference( poly1, poly2, output.surfaces() );
		}
		else {
			output.addPrimitive( surf );
		}
	}

	// prim intersects surf and they share a common plane
	static void substract_from_surface( const PrimitiveBase<3>& surf, const PrimitiveBase<3>& prim, GeometrySet<3>& output )
	{
		if ( prim.getType() == PrimitiveSurface ) {
			const CGAL::Triangle_3<Kernel>& tri1 = surf.as<Surface_d<3>::Type>();
			const CGAL::Triangle_3<Kernel>& tri2 = prim.as<Surface_d<3>::Type>();

			// common plane
			CGAL::Plane_3<Kernel> plane = tri1.supporting_plane();

			CGAL::Polygon_2<Kernel> poly1, poly2;

			poly1.push_back( plane.to_2d( tri1.vertex(0) ) );
			poly1.push_back( plane.to_2d( tri1.vertex(1) ) );
			poly1.push_back( plane.to_2d( tri1.vertex(2) ) );

			poly2.push_back( plane.to_2d( tri2.vertex(0) ) );
			poly2.push_back( plane.to_2d( tri2.vertex(1) ) );
			poly2.push_back( plane.to_2d( tri2.vertex(2) ) );

			PolygonList polys;
			CGAL_poly_difference( poly1, poly2, polys );

			for ( PolygonList::const_iterator it = polys.begin();
			      it != polys.end();
			      ++it ) {
				// add the projected-back polygon
				Polygon p( *it, plane );
				output.addGeometry( p );
			}
		}
		else {
			output.addPrimitive( surf );
		}
	}

	static void substract_from_volume( const PrimitiveBase<2>& vol, const PrimitiveBase<2>& prim, GeometrySet<2>& output )
	{
	}

	static void substract_from_volume( const PrimitiveBase<3>& vol, const PrimitiveBase<3>& prim, GeometrySet<3>& output )
	{
		if ( prim.getType() == PrimitiveVolume ) {
			const MarkedPolyhedron& vol1 = vol.as<Volume_d<3>::Type>();
			const MarkedPolyhedron& vol2 = prim.as<Volume_d<3>::Type>();

			typedef CGAL::Polyhedron_corefinement<MarkedPolyhedron> Corefinement;
			MarkedPolyhedron& polya = const_cast<MarkedPolyhedron&>( vol1 );
			MarkedPolyhedron& polyb = const_cast<MarkedPolyhedron&>( vol2 );

			Corefinement coref;
			CGAL::Emptyset_iterator no_polylines;
			// vector of <Polyhedron, tag>
			std::vector<std::pair<MarkedPolyhedron*, int> > result;
			coref( polya, polyb, no_polylines, std::back_inserter(result), Corefinement::P_minus_Q_tag );

			// empty difference
			if (result.size() == 0) {
				return ;
			}
			
			// else, we have a difference
			MarkedPolyhedron* res_poly = result[0].first;
			output.addPrimitive( *res_poly );
			delete res_poly;
		}
		else {
			output.addPrimitive( vol );
		}
	}

	template <int Dim>
	static void substract_from_primitive( const PrimitiveBase<Dim>& a, const PrimitiveBase<Dim>& b, GeometrySet<Dim>& output )
	{
		switch ( a.getType() ) {
		case PrimitivePoint:
			break;
		case PrimitiveSegment:
			substract_from_segment( a, b, output );
			break;
		case PrimitiveSurface:
			substract_from_surface( a, b, output );
			break;
		case PrimitiveVolume:
			substract_from_volume( a, b, output );
			break;
		}
	}

	template <int Dim>
	// pa and pb intersects
	static void difference_primitive( const PrimitiveBase<Dim>& pa, const PrimitiveBase<Dim>& pb, GeometrySet<Dim>& output )
	{
		if ( pa.getType() == PrimitivePoint ) {
			// difference = empty
		}
		else {
			int paType = pa.getType();

			GeometrySet<Dim> inter;
			algorithm::intersection( pa, pb, inter );

			// FIXME replace by intersectionDim ?
			if ( inter.size( paType ) > 0 ) {
				GeometrySet<Dim> temp( /* isComplete */ true );
				for ( typename GeometrySet<Dim>::const_iterator it = inter.primitives_begin( paType );
				      it != inter.primitives_end();
				      ++it ) {
					GeometrySet<Dim> temp2, temp3;
					substract_from_primitive( pa, *it, temp2 );
					algorithm::intersection( temp, temp2, temp3 );
					temp = temp3;
				}
				output.merge( temp, paType );
			}
			else {
				output.addPrimitive( pa );
			}
		}
	}
	
	template <int Dim>
	void difference( const GeometrySet<Dim>& a, const GeometrySet<Dim>& b, GeometrySet<Dim>& output )
	{
		typename SFCGAL::HandleCollection<Dim>::Type ahandles, bhandles;
		typename SFCGAL::BoxCollection<Dim>::Type aboxes, bboxes;
		a.computeBoundingBoxes( ahandles, aboxes );
		b.computeBoundingBoxes( bhandles, bboxes );


		//
		// (A1 U A2) - B  <=> (A1 - B) U (A2 - B)

		for ( size_t i = 0; i < aboxes.size(); ++i ) {
			// start with the whole plane/volume
			GeometrySet<Dim> temp1( /* isComplete */ true );
			// A - (B1 U B2 ) <=> (A - B1) ^ (A - B2)
			const PrimitiveBase<Dim>* pa = aboxes[i].handle();
			for ( size_t j = 0; j < bboxes.size(); ++j ) {
				GeometrySet<Dim> temp2;
				const PrimitiveBase<Dim>* pb = bboxes[j].handle();
				if ( CGAL::do_overlap(aboxes[i].bbox(), bboxes[j].bbox()) &&
				     algorithm::intersects( *pa, *pb ) ) {
					difference_primitive( *pa, *pb, temp2 );
				}
				else {
					temp2.addPrimitive( *aboxes[i].handle() );
				}

				// merge with temp1
				GeometrySet<Dim> temp3;
				algorithm::intersection( temp1, temp2, temp3 );
				temp1 = temp3;
			}
			// merge with output (only select primitives of A's type)
			output.merge( temp1, pa->getType() );
		}
	}

	template void difference<2>( const GeometrySet<2>& a, const GeometrySet<2>& b, GeometrySet<2>& );
	template void difference<3>( const GeometrySet<3>& a, const GeometrySet<3>& b, GeometrySet<3>& );

	std::auto_ptr<Geometry> difference( const Geometry& ga, const Geometry& gb )
	{
		GeometrySet<2> gsa( ga ), gsb( gb ), output;
		algorithm::difference( gsa, gsb, output );

		return output.recompose();
	}

	std::auto_ptr<Geometry> difference3D( const Geometry& ga, const Geometry& gb )
	{
		GeometrySet<3> gsa( ga ), gsb( gb ), output;
		algorithm::difference( gsa, gsb, output );

		return output.recompose();
	}
}
}
