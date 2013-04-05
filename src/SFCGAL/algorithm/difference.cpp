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
#include <CGAL/Boolean_set_operations_2.h>

#include <SFCGAL/algorithm/difference.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/Geometry.h>
#include <SFCGAL/GeometrySet.h>

namespace SFCGAL {

namespace algorithm
{
	// get the point position on segment
	// in [0,1] if its on the interior
	static Kernel::FT point_position( const CGAL::Segment_2<Kernel>& seg, const CGAL::Point_2<Kernel>& pt )
	{
		CGAL::Point_2<Kernel> pA = seg.source();
		CGAL::Point_2<Kernel> pB = seg.target();
		Kernel::FT num = pt.x() - pA.x();
		Kernel::FT den = pB.x() - pA.x();
		if ( den == 0 ) {
			num = pt.y() - pA.y();
			den = pB.y() - pA.y();
		}
		return num / den;
	}

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

	static void substract_from_surface( const PrimitiveBase<2>& surf, const PrimitiveBase<2>& prim, GeometrySet<2>& output )
	{
		if ( prim.getType() == PrimitiveSurface ) {
			const Surface_d<2>::Type& poly1 = surf.as<Surface_d<2>::Type>();
			const Surface_d<2>::Type& poly2 = prim.as<Surface_d<2>::Type>();

			CGAL::difference( poly1, poly2, std::back_inserter( output.surfaces() ) );
		}
		else {
			output.addPrimitive( surf );
		}
	}

	static void substract_from_surface( const PrimitiveBase<3>& surf, const PrimitiveBase<3>& prim, GeometrySet<3>& output )
	{
#if 0
		if ( prim.getType() == PrimitiveSurface ) {
			Surface_d<3>::Type& tri1 = surf.as<Surface_d<3>::Type>();
			Surface_d<3>::Type& tri2 = prim.as<Surface_d<3>::Type>();

			CGAL::Object interObj = CGAL::intersection( tri1, tri2 );
			const Surface_d<3>::Type* interTri;
			const std::vector<Point_d<3>::Type>* interPts;

			CGAL::Polygon_
			if ( interTri = CGAL::object_cast<Surface_d<3>::Type>( interObj ) ) {
				// the intersection is a triangle,
				
			}
		}
		else {
			output.addPrimitive( surf );
		}
#endif
	}

	template <int Dim>
	// pa and pb intersects
	static void difference_primitive( const PrimitiveBase<Dim>& pa, const PrimitiveBase<Dim>& pb, GeometrySet<Dim>& output )
	{
		if ( pa.getType() == PrimitivePoint ) {
			// difference = empty
		}
		else if ( pa.getType() == PrimitiveSegment ) {
			GeometrySet<Dim> inter;
			algorithm::intersection( pa, pb, inter );
			if ( inter.size( PrimitiveSegment ) > 0 ) {
				for ( typename GeometrySet<Dim>::const_iterator it = inter.primitives_begin( PrimitiveSegment );
				      it != inter.primitives_end();
				      ++it ) {
					substract_from_segment( pa, *it, output );
				}
			}
			else {
				output.addPrimitive( pa );
			}
		}
		else if ( pa.getType() == PrimitiveSurface ) {
			substract_from_surface( pa, pb, output );
		}
	}
	
	template <int Dim>
	static void filter_self_intersection( const GeometrySet<Dim>& input, GeometrySet<Dim>& output )
	{
		for ( int ptype = PrimitivePoint; ptype < PrimitiveVolume; ++ptype ) {
			GeometrySet<Dim> temp;

			// idx is a table of pointer to primitives
			std::vector< const PrimitiveBase<Dim>* > idx;

			// first consider input primitives
			for ( typename GeometrySet<Dim>::const_iterator it = input.primitives_begin( ptype );
			      it != input.primitives_end();
			      ++it ) {
				idx.push_back( &*it );
			}

			for ( size_t i = 0; i < idx.size(); ++i ) {
				if ( !idx[i] ) continue;
				bool intersectsA = false;
				for ( size_t j = i+1; j < idx.size(); ++j ) {
					if ( !idx[j] ) continue;
					if ( CGAL::do_overlap( idx[i]->bbox(), idx[j]->bbox() ) &&
					     algorithm::intersects( *idx[i], *idx[j] ) ) {
						intersectsA = true;
						algorithm::intersection( *idx[i], *idx[j], temp );

						for ( typename GeometrySet<Dim>::const_iterator it = temp.primitives_begin( ptype );
						      it != temp.primitives_end();
						      ++it ) {
							// add intersection primitives to the table
							idx.push_back( &*it );
						}
						// the jth primitive has been replaced by an intersection,
						// remove it from the table
						idx[j] = 0;
						break;
					}
				}
				if ( ! intersectsA ) {
					output.addPrimitive( *idx[i] );
				}
				idx[i] = 0;
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

		for ( size_t i = 0; i < aboxes.size(); ++i ) {
			bool intersectsA = false;
			GeometrySet<Dim> tempOut;
			for ( size_t j = 0; j < bboxes.size(); ++j ) {
				if ( CGAL::do_overlap(aboxes[i].bbox(), bboxes[j].bbox()) ) {
					const PrimitiveBase<Dim>* pa = aboxes[i].handle();
					const PrimitiveBase<Dim>* pb = bboxes[j].handle();
					if ( algorithm::intersects( *pa, *pb ) ) {
						intersectsA = true;

						difference_primitive( *pa, *pb, tempOut );
					}
				}
			}
			if ( ! intersectsA ) {
				tempOut.addPrimitive( *aboxes[i].handle() );
			}

			filter_self_intersection( tempOut, output );
		}
	}

	template void difference<2>( const GeometrySet<2>& a, const GeometrySet<2>& b, GeometrySet<2>& );
	//	template void difference<3>( const GeometrySet<3>& a, const GeometrySet<3>& b, GeometrySet<3>& );

	std::auto_ptr<Geometry> difference( const Geometry& ga, const Geometry& gb )
	{
		GeometrySet<2> gsa( ga ), gsb( gb ), output;
		algorithm::difference( gsa, gsb, output );

		return output.recompose();
	}

}
}
