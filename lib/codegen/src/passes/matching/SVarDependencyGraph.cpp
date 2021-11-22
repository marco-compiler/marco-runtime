#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/detail/adjacency_list.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <cstddef>
#include <llvm/ADT/ArrayRef.h>
#include <llvm/ADT/Optional.h>
#include <llvm/ADT/SmallVector.h>
#include <llvm/InitializePasses.h>
#include <map>
#include <marco/mlirlowerer/passes/matching/MatchedEquationLookup.h>
#include <marco/mlirlowerer/passes/matching/SVarDependencyGraph.h>
#include <marco/mlirlowerer/passes/matching/VVarDependencyGraph.h>
#include <marco/mlirlowerer/passes/model/Model.h>
#include <marco/mlirlowerer/passes/model/Variable.h>
#include <optional>

using namespace marco::codegen::model;
using namespace llvm;
using namespace boost;
using namespace std;

void SVarDepencyGraph::insertNode(LookUp& LookUp, size_t vertexIndex)
{
	const IndexesOfEquation& vertex = collapsedGraph[vertexIndex];
	const auto& interval = vertex.getEquation().getInductions();

	for (auto eqInds : interval.contentRange())
	{
		auto indicies = vertex.getEqToVar().map(eqInds);
		auto vertexIndex = add_vertex(SingleEquationReference(vertex, indicies), graph);
		LookUp[&vertex][vertex.getVariable().indexOfElement(indicies)] = vertexIndex;
	}
}

static Optional<size_t> indexOfScalarVar(
		ArrayRef<size_t> access,
		const IndexesOfEquation& var,
		const SVarDepencyGraph::LookUp& lookUp)
{
	auto v = lookUp.find(&var);

	if (v == lookUp.end())
		return {};

	auto toReturn = v->second.find(var.getVariable().indexOfElement(access));

	if (toReturn == v->second.end())
		return {};

	return toReturn->second;
}

void SVarDepencyGraph::insertEdge(
		const LookUp& lookUp, const VVarDependencyGraph::EdgeDesc& edge)
{
	size_t sourceVertex = source(edge, collImpl());
	size_t targetVertex = target(edge, collImpl());
	const IndexesOfEquation& targetNode = *collImpl()[targetVertex];

	if (lookUp.find(&targetNode) == lookUp.end())
		return;

	const IndexesOfEquation& sourceNode = *collImpl()[sourceVertex];
	const VectorAccess& varAccess = collImpl()[edge];

	VectorAccess dependencies = targetNode.getVarToEq() * varAccess;

	for (const auto& indecies : targetNode.getInterval().contentRange())
	{
		auto sourceIndex = indexOfScalarVar(indecies, sourceNode, lookUp);

		auto scalarVarInduction = dependencies.map(indecies);
		auto targetIndex = indexOfScalarVar(scalarVarInduction, targetNode, lookUp);

		if (!sourceIndex || !targetIndex)
			continue;

		add_edge(*sourceIndex, *targetIndex, graph);
	}
}

void SVarDepencyGraph::insertEdges(const LookUp& lookUp, size_t vertexIndex)
{
	for (const auto& edge : outEdgesRange(vertexIndex, collImpl()))
		insertEdge(lookUp, edge);
}

SVarDepencyGraph::SVarDepencyGraph(
		const VVarDependencyGraph& collapsedGraph, const VVarScc& scc)
		: scc(scc), collapsedGraph(collapsedGraph)
{
	LookUp vertexesLookUp;

	for (const auto& vertex : scc)
		insertNode(vertexesLookUp, vertex);

	for (size_t vertex : scc)
		insertEdges(vertexesLookUp, vertex);
}