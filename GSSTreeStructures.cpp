#include "GSSTreeStructures.h"
#include "DataViewers.h"

//write out a leaf node to outNodes and the representative trees to outTrees
void GSSLeaf::writeLeaf(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees)
{
	GSSNodeCommon info;
	info.isLeaf = true;
	info.N = cluster.size();

	unsigned positions[info.N];

	outNodes.write((char*)&info, sizeof(info));

	//compute the position offsets
	unsigned curoffset = info.N*sizeof(unsigned); //offsets are from start of data section
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		positions[i] = curoffset;
		const MappableOctTree *tree = data->getMIV(ind);
		curoffset += sizeof(file_index);
		curoffset += tree->bytes();
	}

	//write positions
	outNodes.write((char*)positions, sizeof(unsigned)*info.N);

	//write tree info
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		//for leaf, write out object file index
		file_index index = data->getIndex(ind);
		outNodes.write((char*)&index, sizeof(file_index));

		//and single tree
		const MappableOctTree *tree = data->getMIV(ind);
		tree->write(outNodes);
	}

	//write out trees of cluster to a GSSDoubleTree
	unsigned offset = cluster.MIV->bytes();

	outTrees.write((char*)&offset, sizeof(offset));
	cluster.MIV->write(outTrees);
	cluster.MSV->write(outTrees);

}

//very similar to write leaf, but write out MIV/MSV of subnodes
void GSSInternalNode::writeNode(const DataViewer *data, const Cluster& cluster, ostream& outNodes, ostream& outTrees)
{
	GSSNodeCommon info;
	info.isLeaf = true;
	info.N = cluster.size();

	unsigned positions[info.N];

	outNodes.write((char*)&info, sizeof(info));

	//compute the position offsets
	unsigned curoffset = info.N*sizeof(unsigned); //offsets are from start of data section
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		positions[i] = curoffset;
		const MappableOctTree *miv = data->getMIV(ind);
		const MappableOctTree *msv = data->getMIV(ind);
		curoffset += sizeof(file_index);
		curoffset += sizeof(unsigned);
		curoffset += miv->bytes();
		curoffset += msv->bytes();
	}

	//write positions
	outNodes.write((char*)positions, sizeof(unsigned)*info.N);

	//write tree info
	for(unsigned i = 0, n = cluster.size(); i < n; i++)
	{
		unsigned ind = cluster[i];
		//for leaf, write out node file index
		file_index index = data->getIndex(ind);
		outNodes.write((char*)&index, sizeof(file_index));

		//and both trees
		const MappableOctTree *miv = data->getMIV(ind);
		const MappableOctTree *msv = data->getMIV(ind);
		unsigned offset = miv->bytes();
		outNodes.write((char*)&offset, sizeof(offset));
		miv->write(outNodes);
		msv->write(outNodes);
	}

	//write out trees of cluster to a GSSDoubleTree
	unsigned offset = cluster.MIV->bytes();

	outTrees.write((char*)&offset, sizeof(offset));
	cluster.MIV->write(outTrees);
	cluster.MSV->write(outTrees);

}
