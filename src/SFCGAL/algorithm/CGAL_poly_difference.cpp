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

#include <SFCGAL/Kernel.h>
#include <SFCGAL/GeometrySet.h>

#include <CGAL/Boolean_set_operations_2.h>

namespace SFCGAL
{

typedef std::list<CGAL::Polygon_with_holes_2<Kernel> > PolygonList;

void CGAL_poly_difference( const CGAL::Polygon_2<Kernel> p1, const CGAL::Polygon_2<Kernel> p2, PolygonList& l)
{
	CGAL::difference( p1, p2, std::back_inserter( l ) );
}

void CGAL_poly_difference( const CGAL::Polygon_with_holes_2<Kernel> p1, const CGAL::Polygon_with_holes_2<Kernel> p2, GeometrySet<2>::SurfaceCollection& l)
{
	CGAL::difference( p1, p2, std::back_inserter( l ) );
}

}
