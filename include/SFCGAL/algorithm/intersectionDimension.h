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
#ifndef SFCGAL_INTERSECTION_DIM_ALGORITHM
#define SFCGAL_INTERSECTION_DIM_ALGORITHM

namespace SFCGAL {
	class Geometry;
	template <int Dim> class GeometrySet;
	template <int Dim> class PrimitiveBase;

	namespace algorithm {
	/*
	 * Intersection dimension between two 2D geometries
	 * -1: no intersection
	 *  0: the intersection is a point
	 *  1: the intersection is a segment/line
	 *  2: the intersection is a surface
	 */
	int intersectionDimension( const Geometry& ga, const Geometry& gb );

	template <int Dim>
	int intersectionDimension( const GeometrySet<Dim>& a, const GeometrySet<Dim>& b );

	template <int Dim>
	int intersectionDimension( const PrimitiveBase<Dim>& a, const PrimitiveBase<Dim>& b );

    }
}

#endif
