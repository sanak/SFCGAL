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

#ifndef _SFCGAL_DETAIL_GEOMETRY_SET_H_
#define _SFCGAL_DETAIL_GEOMETRY_SET_H_

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/variant.hpp>

#include <SFCGAL/Kernel.h>
#include <SFCGAL/TypeForDimension.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Box_intersection_d/Box_with_handle_d.h>

// comparison operator on segments, for use in a std::set
bool operator< ( const CGAL::Segment_2<SFCGAL::Kernel>& sega, const CGAL::Segment_2<SFCGAL::Kernel>& segb );
bool operator< ( const CGAL::Segment_3<SFCGAL::Kernel>& sega, const CGAL::Segment_3<SFCGAL::Kernel>& segb );

namespace SFCGAL {
	class Geometry;

	///
	/// Primitive type enumeration.
	/// The associated integer is the dimension of the primitive
	enum PrimitiveType
	{
		PrimitivePoint,
		PrimitiveSegment,
		PrimitiveSurface,
		PrimitiveVolume
	};

	///
	/// Flags available for each type of Primitive.
	/// Primitives can be 'flagged' in order to speed up recomposition
	enum PrimitiveFlag
	{
		// the polyhedron is planar => build a triangle or a polygon
		FLAG_IS_PLANAR = 1
	};

	template <int Dim, class T> class Primitive;
	///
	/// Base class for Primitive. Stores a flag
	template <int Dim>
	class PrimitiveBase
	{
	public:
		PrimitiveBase() : _flags(0) {}
		PrimitiveBase( int f ) : _flags(f) {}

		///
		/// get the primitive type ( point, segment, surface, volume )
		virtual int getType() const = 0;

		int flags() const { return _flags; }
		void setFlags( int flags ) { _flags = flags; }

		///
		/// Access to the primitive
		/// Example usage: as<Point_d<2>::Type>()
		template <class T>
		const T& as() const {
			return static_cast<const Primitive<Dim, T>* >( this )->primitive();
		}
		template <class T>
		T& as() {
			return static_cast<Primitive<Dim, T>* >( this )->primitive();
                }

	        ///
	        /// Compute the bounding box
		virtual typename Bbox_d<Dim>::Type bbox() const = 0;
	private:
		int _flags;
	};

	///
	/// Primitive class
	/// T : Point_d, Segment_d, Surface_d, Volume_d
        template < int Dim, class T >
	class Primitive : public PrimitiveBase<Dim>
	{
	public:
		Primitive() :  PrimitiveBase<Dim>() {}
		Primitive( const T& p ) : PrimitiveBase<Dim>(), _primitive(p) {}
		Primitive( const T& p, int f ) : PrimitiveBase<Dim>(f), _primitive(p) {}

		bool operator< ( const Primitive& other ) const {
			return _primitive < other._primitive;
		}

		int getType() const {
			return PrimitiveDimension<T>::value;
		}

		T& primitive() { return _primitive; }
		const T& primitive() const { return _primitive; }

		typename Bbox_d<Dim>::Type bbox() const {
			return computeBbox( *this );
		}
	private:
		T _primitive;
	};

	template <int Dim>
	struct PrimitivePoint_d
	{
		typedef Primitive<Dim, typename Point_d<Dim>::Type> Type;
	};
	template <int Dim>
	struct PrimitiveSegment_d
	{
		typedef Primitive<Dim, typename Segment_d<Dim>::Type> Type;
	};
	template <int Dim>
	struct PrimitiveSurface_d
	{
		typedef Primitive<Dim, typename Surface_d<Dim>::Type> Type;
	};
	template <int Dim>
	struct PrimitiveVolume_d
	{
		typedef Primitive<Dim, typename Volume_d<Dim>::Type> Type;
	};

        template <int Dim, class T>
	std::ostream& operator<<( std::ostream& ostr, const Primitive<Dim, T>& p )
	{
		ostr << p.primitive() << " F:" << p.flags();
		return ostr;
	}

        template <int Dim, class T>
        typename Bbox_d<Dim>::Type computeBbox( const Primitive<Dim, T>& p )
        {
		return p.primitive().bbox();
        }

        Bbox_d<2>::Type computeBbox( const PrimitiveVolume_d<2>::Type& p );
        Bbox_d<3>::Type computeBbox( const PrimitiveVolume_d<3>::Type& p );

	///
	/// PrimitiveBox. Type used for CGAL::Box_intersection_d
	template <int Dim>
	struct PrimitiveBox
	{
		typedef CGAL::Box_intersection_d::Box_with_handle_d<double, Dim, const PrimitiveBase<Dim>*> Type;
	};

	///
	/// BoxCollection for use with CGAL::Box_intersection_d
	template <int Dim>
	struct BoxCollection
	{
		typedef std::vector<typename PrimitiveBox<Dim>::Type> Type;
	};

        ///
	/// HandleCollection. Used to store PrimitiveHandle
        template <int Dim>
        struct HandleCollection
	{
		typedef std::list<const PrimitiveBase<Dim>* > Type;
	};

        ///
        /// Used to choose between a regular iterator and a const_iterator (needed for GeometrySet<>::Iterator
	template <class T, bool IsConst>
	struct IteratorChooser {
		typedef typename T::const_iterator Type;
	};
	template <class T>
	struct IteratorChooser<T, false> {
		typedef typename T::iterator Type;
	};

	///
	/// A GeometrySet represents a set of CGAL primitives.
	/// Primitive are either of dimension 0 (points),
	/// dimension 1 (segments), dimension 2 (surfaces, a.k.a. polygon or triangles)
	/// or dimension 3 (polyhedron)
	template <int Dim>
	class GeometrySet
	{
	public:
		// Points are stored in an ordered set
		typedef std::set<typename PrimitivePoint_d<Dim>::Type > PointCollection;
		// Segments are stored in an ordered set
		typedef std::set<typename PrimitiveSegment_d<Dim>::Type > SegmentCollection;
		typedef std::list<typename PrimitiveSurface_d<Dim>::Type > SurfaceCollection;
		typedef std::list<typename PrimitiveVolume_d<Dim>::Type > VolumeCollection;

		template <class T>
		class IteratorBase : public boost::iterator_facade < IteratorBase<T>,
								     T,
								     boost::forward_traversal_tag >
		{
		public:
			typedef typename IteratorChooser< PointCollection, boost::is_const<T>::value >::Type PointIterator;
			typedef typename IteratorChooser< SegmentCollection, boost::is_const<T>::value >::Type SegmentIterator;
			typedef typename IteratorChooser< SurfaceCollection, boost::is_const<T>::value >::Type SurfaceIterator;
			typedef typename IteratorChooser< VolumeCollection, boost::is_const<T>::value >::Type VolumeIterator;

			IteratorBase() {}
			IteratorBase( PointIterator pit,
				      PointIterator pend,
				      SegmentIterator lit,
				      SegmentIterator lend,
				      SurfaceIterator sit,
				      SurfaceIterator send,
				      VolumeIterator vit,
				      VolumeIterator vend );
			
			void increment();

			bool equal( const IteratorBase<T>& other ) const;

			T& dereference() const;

		private:
			PointIterator pointIt, pointEnd;
			SegmentIterator segmentIt, segmentEnd;
			SurfaceIterator surfaceIt, surfaceEnd;
			VolumeIterator volumeIt, volumeEnd;

			friend IteratorBase<T> GeometrySet<Dim>::insert( IteratorBase<T>, const T& );
		};

		typedef IteratorBase<const PrimitiveBase<Dim> > ConstIterator;
		typedef IteratorBase<PrimitiveBase<Dim> >       Iterator;

		GeometrySet();

		/**
		 * Construct a GeometrySet from a SFCGAL::Geometry
		 */
		GeometrySet( const Geometry& g );

		/**
		 * Construct a GeometrySet from a Point
		 */
		GeometrySet( const typename PrimitivePoint_d<Dim>::Type& g );

		/**
		 * Construct a GeometrySet from a Segment
		 */
		GeometrySet( const typename PrimitiveSegment_d<Dim>::Type& g );

		/**
		 * Construct a GeometrySet from a Surface
		 */
		GeometrySet( const typename PrimitiveSurface_d<Dim>::Type& g );

		/**
		 * Construct a GeometrySet from a Volume
		 */
		GeometrySet( const typename PrimitiveVolume_d<Dim>::Type& g );

		/**
		 * Add a geometry by decomposing it into CGAL primitives
		 */
		void addGeometry( const Geometry& g );

		/**
		 * add a primitive from a PrimitiveHandle  to the set
		 */
		void addPrimitive( const PrimitiveBase<Dim>& p );

		/**
		 * add a primitive from a CGAL::Object to the set
		 * pointsAsRing : if set to true, build a polygon if o is a vector of points
		 */
		void addPrimitive( const CGAL::Object& o, bool pointsAsRing = false );

		/**
		 * add a point to the set
		 */
		void addPrimitive( const typename PrimitivePoint_d<Dim>::Type& g );
		template <class IT>
		void addPoints( IT ibegin, IT iend )
		{
			std::copy( ibegin, iend, std::inserter(_points, _points.end()) );
		}

		/**
		 * collect all points of b and add them to the point list
		 */
		void collectPoints( const PrimitiveBase<Dim>& b );

		/**
		 * add a segment to the set
		 */
		void addPrimitive( const typename PrimitiveSegment_d<Dim>::Type& g );
		template <class IT>
		void addSegments( IT ibegin, IT iend )
		{
			std::copy( ibegin, iend, std::inserter(_segments, _segments.end()) );
		}

		/**
		 * add a surface to the set
		 */
		void addPrimitive( const typename PrimitiveSurface_d<Dim>::Type& g );
		template <class IT>
		void addSurfaces( IT ibegin, IT iend )
		{
			std::copy( ibegin, iend, std::back_inserter(_surfaces) );
		}

		/**
		 * add a volume to the set
		 */
		void addPrimitive( const typename PrimitiveVolume_d<Dim>::Type& g );
		template <class IT>
		void addVolumes( IT ibegin, IT iend )
		{
			std::copy( ibegin, iend, std::back_inserter(_volumes) );
		}

		/**
		 * Compute all bounding boxes and handles of the set
		 */
		void computeBoundingBoxes( typename HandleCollection<Dim>::Type& handles, typename BoxCollection<Dim>::Type& boxes ) const;

		inline PointCollection& points() { return _points; }
		inline const PointCollection& points() const { return _points; }

		inline SegmentCollection& segments() { return _segments; }
		inline const SegmentCollection& segments() const { return _segments; }

		inline SurfaceCollection& surfaces() { return _surfaces; }
		inline const SurfaceCollection& surfaces() const { return _surfaces; }

		inline VolumeCollection& volumes() { return _volumes; }
		inline const VolumeCollection& volumes() const { return _volumes; }

		///
		/// Primitive iterator
		/// the selection parameter is set to the proper PrimitiveType the caller wants to iterate over
		/// -1 is passed to iterate over all the different kinds of primitives
		ConstIterator primitives_begin( int selection = -1 ) const;
		ConstIterator primitives_end() const;

		Iterator primitives_begin( int selection = -1 );
		Iterator primitives_end();

		/**
		 * STL-compliant insert
		 * For use with std::inserter( gset, it )
		 */
		Iterator insert( Iterator it, const PrimitiveBase<Dim>& value );
		ConstIterator insert( ConstIterator it, const PrimitiveBase<Dim>& value );

		/**
		 * convert the set to a SFCGAL::Geometry
		 */
		std::auto_ptr<Geometry> recompose() const;

		/**
		 * Filter (remove) primitives that are already covered by others
		 */
		// TODO : move out of the class
		void filterCovered( GeometrySet<Dim>& output ) const;

	private:
		///
		/// Given an input SFCGAL::Geometry, decompose it into CGAL primitives
		void _decompose( const Geometry& g );

		PointCollection _points;
		SegmentCollection _segments;
		SurfaceCollection _surfaces;
		VolumeCollection _volumes;
	};

	///
	/// Display operator
	std::ostream& operator<<( std::ostream&, const GeometrySet<2>& g );
	///
	/// Display operator
	std::ostream& operator<<( std::ostream&, const GeometrySet<3>& g );
}

#endif
