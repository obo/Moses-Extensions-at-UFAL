#include "ChartTrellisPath.h"
#include "ChartHypothesis.h"
#include "ChartTrellisPathCollection.h"
#include "../../moses/src/StaticData.h"
#include "../../moses/src/Word.h"

using namespace std;

namespace MosesChart
{

TrellisPath::TrellisPath(const Hypothesis *hypo)
:m_finalNode(new TrellisNode(hypo))
,m_scoreBreakdown(hypo->GetScoreBreakDown())
,m_totalScore(hypo->GetTotalScore())
,m_prevNodeChanged(NULL)
,m_prevPath(NULL)
{
}

TrellisPath::TrellisPath(const TrellisPath &origPath
												, const TrellisNode &soughtNode
												, const Hypothesis &replacementHypo
												, Moses::ScoreComponentCollection	&scoreChange)
:m_scoreBreakdown(origPath.GetScoreBreakdown())
,m_prevPath(&origPath)
{
	m_finalNode = new TrellisNode(origPath.GetFinalNode()
															, soughtNode
															, replacementHypo
															, scoreChange
															, m_prevNodeChanged);

	m_scoreBreakdown.PlusEquals(scoreChange);

	m_totalScore = m_scoreBreakdown.GetWeightedScore();
}

TrellisPath::~TrellisPath()
{
	delete m_finalNode;
}

Moses::Phrase TrellisPath::GetOutputPhrase() const
{
	Moses::Phrase ret = GetFinalNode().GetOutputPhrase();
	return ret;
}

void TrellisPath::CreateDeviantPaths(TrellisPathCollection &pathColl, const TrellisNode &soughtNode) const
{
	// copy this path but replace startHypo with its arc
	const ArcList *arcList = soughtNode.GetHypothesis().GetArcList();

	if (arcList)
	{
		ArcList::const_iterator iterArcList;
		for (iterArcList = arcList->begin(); iterArcList != arcList->end(); ++iterArcList)
		{
			Moses::ScoreComponentCollection	scoreChange;

			const Hypothesis &replacementHypo = **iterArcList;
			TrellisPath *newPath = new TrellisPath(*this, soughtNode, replacementHypo, scoreChange);
			pathColl.Add(newPath);
		}
	}

	// recusively create deviant paths for child nodes
	const TrellisNode::NodeChildren &children = soughtNode.GetChildren();

	TrellisNode::NodeChildren::const_iterator iter;
	for (iter = children.begin(); iter != children.end(); ++iter)
	{
		const TrellisNode &child = **iter;
		CreateDeviantPaths(pathColl, child);
	}
}

void TrellisPath::CreateDeviantPaths(TrellisPathCollection &pathColl) const
{
	if (m_prevNodeChanged == NULL)
	{ // initial enumeration from a pure hypo
		CreateDeviantPaths(pathColl, GetFinalNode());
	}
	else
	{ // don't change m_prevNodeChanged, just it's children
		const TrellisNode::NodeChildren &children = m_prevNodeChanged->GetChildren();

		TrellisNode::NodeChildren::const_iterator iter;
		for (iter = children.begin(); iter != children.end(); ++iter)
		{
			const TrellisNode &child = **iter;
			CreateDeviantPaths(pathColl, child);
		}
	}
}

std::ostream& operator<<(std::ostream &out, const TrellisPath &path)
{
	out << &path << "  " << path.m_prevPath << "  " << path.GetOutputPhrase() << endl;

	if (path.m_prevNodeChanged)
	{
		out << "changed " << path.m_prevNodeChanged->GetHypothesis() << endl;
	}

	out << path.GetFinalNode() << endl;

	return out;
}

};

