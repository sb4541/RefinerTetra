#ifndef __CUTTER_3D_H
#define __CUTTER_3D_H

#include <vector>
#include <string>
#include <set>
#include <list>

#include "GenericDomain.hpp"
#include "GenericMesh.hpp"
#include "Intersector2D1D.hpp"
#include "Output.hpp"
#include "Eigen/Eigen"

using namespace std;
using namespace MainApplication;

namespace GeDiM
{
    class CutterByDomains;
	class CutterMesh3D;

    class CutterByDomains
	{
		protected:
			GenericMesh* mesh;
			vector<GenericDomain*> domains;
			map<unsigned int, unsigned int> positionDomains;

		public:
			CutterByDomains();
			virtual ~CutterByDomains();

			const Output::ExitCodes AddDomain(GenericDomain& domain);
			const Output::ExitCodes Initialize(const int& numberOfDomain);

			const Output::ExitCodes SetMesh(GenericMesh& _mesh) {mesh = &_mesh; return Output::Success;}
			virtual const Output::ExitCodes CutMesh() = 0;

			const GenericDomain& DomainByPosition(const unsigned int& position) const {return *domains[position];}
			GenericDomain& DomainByPosition(const unsigned int& position) {return *domains[position];}
			const GenericDomain& DomainById(const unsigned int& id) const {return *domains[positionDomains.at(id)];}
			GenericDomain& DomainById(const unsigned int& id) {return *domains[positionDomains.at(id)];}

			void SortDomainsBySegmentOnBorder(const GenericDomain2D* domain);
			void SortDomainsByMeasure(const bool& ascendent = true){if(!ascendent) sort(domains.begin(),domains.end(), CompareByMeasureDescendent); else sort(domains.begin(),domains.end(), CompareByMeasureAscendent);}
			void SortDomainsByPlane(){sort(domains.begin(),domains.end(), CompareByPlane);}

			static bool CompareByMeasureDescendent(const GenericDomain* first, const GenericDomain* second) { return first->Measure() > second->Measure(); }
			static bool CompareByMeasureAscendent(const GenericDomain* first, const GenericDomain* second) { return first->Measure() < second->Measure(); }
			static bool CompareByPlane(const GenericDomain* first, const GenericDomain* second);
	};

	class CutterMesh3D : public CutterByDomains
	{
		protected:
			enum PositionInRelationToSurface
			{
				Plane = 0,
				Positive = 1,
				Negative = 2,
				ToBreak = 3
			};

			Intersector2D1D* intersector2D1D;

			vector<PositionInRelationToSurface> valuePoints;
			map<unsigned int, unsigned int> idPositionPoints;

			vector<PositionInRelationToSurface> valueEdges;
			map<unsigned int, unsigned int> idPositionEdges;

			set<unsigned int> positivePoints;
			set<unsigned int> negativePoints;
			set<unsigned int> planePoints;

			set<unsigned int> positiveEdges;
			set<unsigned int> negativeEdges;
			set<unsigned int> planeEdges;

			list<unsigned int> positiveFaces;
			list<unsigned int> negativeFaces;
			list<unsigned int> planeFaces;

			set<unsigned int > idFacesToUpdate;
			set<unsigned int > idCellsToUpdate;

			static bool CompareBySizeList(const list<GenericCell*>& first, const list<GenericCell*>& second) {return first.size() > second.size();}

			//metodi di frattura
			const Output::ExitCodes CutFace(GenericFace& face, const vector<const GenericPoint*>& domainPointsInFace, const double& toll = 1.0E-7);
			const Output::ExitCodes CutEdge(GenericEdge& edge, unsigned int& idNewPoint);

			//metodi di controllo rispetto il dominio
			const PositionInRelationToSurface PointDomainSplitter(const GenericPoint& point, const Vector3d& normalPlane, const double& translation, const double& toll = 1.0E-7);
			const PositionInRelationToSurface EdgeDomainSplitter(const GenericEdge& edge);
			const bool FaceDomainSplitter(GenericFace& face, vector<const GenericPoint*>& domainPointsInFace);
			const bool FaceDomainSplitter2(GenericFace& face, vector<const GenericPoint*>& domainPointsInFace);


		public:
			CutterMesh3D();
			virtual ~CutterMesh3D();

			void SetToleranceParallelism(const double& toll) {intersector2D1D->SetToleranceParallelism(toll);}
			void SetToleranceIntersection(const double& toll) {intersector2D1D->SetToleranceIntersection(toll);}
			const double& ToleranceParellism() const { return intersector2D1D->ToleranceParallelism();}
			const double& ToleranceIntersection() const { return intersector2D1D->ToleranceIntersection();}

			const Output::ExitCodes CutCell(GenericCell& cell, const Vector3d& normalPlane, const double& translation, const int& idDomain = -1, const double& toll = 1.0E-7);
			const Output::ExitCodes Reset();

			const Output::ExitCodes CutMesh();
	};
}

#endif // __CUTTER_3D_H

