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

#include <SFCGAL/algorithm/intersectionDimension.h>
#include <SFCGAL/algorithm/intersects.h>
#include <SFCGAL/algorithm/intersection.h>
#include <SFCGAL/GeometrySet.h>

namespace SFCGAL {
namespace algorithm
{
	template <int Dim>
	int intersectionDimension( const PrimitiveBase<Dim>& pa, const PrimitiveBase<Dim>& pb )
	{
		if ( pa.getType() < pb.getType() ) {
			// swap arguments
			return intersectionDimension( pb, pa );
		}

		if ( ! algorithm::intersects( pa, pb ) ) {
			return -1;
		}
		// we have pa.getType() >= pb.getType()
		if ( pb.getType() == PrimitivePoint ) {
			return 0;
		}
		
		if ( pa.getType() == PrimitiveSegment ) {
			const typename Segment_d<Dim>::Type& sega = pa.template as<typename Segment_d<Dim>::Type>();
			if ( pb.getType() == PrimitiveSegment ) {
				const typename Segment_d<Dim>::Type& segb = pb.template as<typename Segment_d<Dim>::Type>();
				typename PrimitivePoint_d<Dim>::Type pas( sega.source() );
				typename PrimitivePoint_d<Dim>::Type pat( sega.target() );
				typename PrimitivePoint_d<Dim>::Type pbs( segb.source() );
				typename PrimitivePoint_d<Dim>::Type pbt( segb.target() );
				if ( algorithm::intersects( pa, pbs ) ||
				     algorithm::intersects( pa, pbt ) ||
				     algorithm::intersects( pb, pas ) ||
				     algorithm::intersects( pb, pat ) ) {
					return 1;
				}
				return 0;
			}
		}
		// default
		GeometrySet<Dim> inter;
		algorithm::intersection( pa, pb, inter );
		return inter.maximumDimension();
	}

	template <int Dim>
	struct intersectiondim_cb
	{
		intersectiondim_cb() : maxdim(-1) {}
		void operator()( const typename PrimitiveBox<Dim>::Type& a,
				 const typename PrimitiveBox<Dim>::Type& b )
		{
			int dim = intersectionDimension( *a.handle(), *b.handle() );
			if ( dim > maxdim ) {
				maxdim = dim;
			}
		}
		int maxdim;
	};

	template <int Dim>
	int intersectionDimension( const GeometrySet<Dim>& a, const GeometrySet<Dim>& b )
	{
		if ( a.complete() ) {
			return b.maximumDimension();
		}
		if ( b.complete() ) {
			return a.maximumDimension();
		}
		typename SFCGAL::HandleCollection<Dim>::Type ahandles, bhandles;
		typename SFCGAL::BoxCollection<Dim>::Type aboxes, bboxes;
		a.computeBoundingBoxes( ahandles, aboxes );
		b.computeBoundingBoxes( bhandles, bboxes );

		intersectiondim_cb<Dim> cb;
		CGAL::box_intersection_d( aboxes.begin(), aboxes.end(),
					  bboxes.begin(), bboxes.end(),
					  cb );
		return cb.maxdim;
	}

	template int intersectionDimension<2>( const PrimitiveBase<2>& pa, const PrimitiveBase<2>& pb );
	template int intersectionDimension<3>( const PrimitiveBase<3>& pa, const PrimitiveBase<3>& pb );

	int intersectionDimension( const Geometry& ga, const Geometry& gb )
	{
		GeometrySet<2> gsa( ga );
		GeometrySet<2> gsb( gb );

		return intersectionDimension( gsa, gsb );
	}
}
}
