#include "SCTree.h"
#include <iostream>

SCTree::SCTree()
{
	// Create a root node for the tree
	root = new SCNode();
	root->name = "root";
	root->parent = NULL;
	nodelist.push_back(root);
}


SCTree::SCTree(string name)
{
	// Create a root node for the tree
	root = new SCNode();
	root->name = name;
	root->parent = NULL;
	nodelist.push_back(root);
}


SCTree::~SCTree()
{
	// Delete Nodes in the Tree
//	for (unsigned int i = 0; i < nodelist.size(); ++i) {
//		if (nodelist[i] != NULL) {		
//			delete nodelist[i];
//			nodelist[i] = NULL;
//		}
//	}
}


//! Rapid display on stdout of the tree structure
void 
SCTree::DrawOnTerminal(bool distances)
{
	cout << "\n\nTerminal Representation of the Tree:\n";
	root->DrawOnTerminal(0, distances);
}


unsigned int
SCTree::FindParent2(string name)
{
 	return parentlist2[name];
}


string 
SCTree::GetTreeString(
	bool distances, 
	double scaleFactor)
{
	string tree;

	GetTreeRecurse(tree, root, distances, scaleFactor);

	return tree + ";";
}


void 
SCTree::GetTreeRecurse(
	string& ret, 
	SCNode* node,
	bool distances, 
	double scaleFactor)
{
	string distance = "";

	if (distances) {
		double bl = node->GetDistance();

		// We don't want to output exponential notation
		//cout << "bl is: " << bl << endl;
		//if (bl < 0.0001)
		//  distance = ":0.0001";
		//else {
		//if (bl > 9)
		//bl = 0.9;
		char buf[1024];
		sprintf(buf, ":%f", bl);
		distance = buf;
			//}
	}

	if (node->IsLeaf()) {
	  //cout << "I am a leaf node!" << endl;
	  //cout << "my name is: " << node-> name << endl;
//		ret.append(node->name + ":" + distance);
		ret.append(node->name);
		return;
	}


	
// CHECK THIS

	unsigned int numChildren = node->NumChildren();
	if (numChildren!= 1)
	  ret += "(";
	if (numChildren == 2) {
	  
	  if (node->children[0]->name > node->children[1]->name) {
	    SCNode *temp = node->children[0];
	    node->children[0] = node->children[1];
	    node->children[1] = temp;
	  }
	  
	  if (!node->children[0]->IsLeaf() && !node->children[1]->IsLeaf()) {
	    SCNode *temp, *temp2;
	    temp = GetLeastSubtree(node);
	    temp2 = temp->parent;
	    
	    while (temp2 != node) {
	      temp = temp2;
	      temp2 = temp2->parent;
	    }
	    
	    if (temp->name == node->children[1]->name) {
	      node->children[1] = node->children[0];
	      node->children[0] = temp;
	    }
	  }
	  
	}
    
	for (unsigned int i = 0; i < numChildren; i++) {

	  GetTreeRecurse(ret, node->children[i], distances, scaleFactor);
	  //		GetTreeRecurse(ret, node->children[i]);
	  if (i != numChildren - 1)
	    ret += ",";

	}
	
	if (node->IsRoot()) {
		ret += ")";
		return;
	}

	if (numChildren!= 1)
	  ret = ret + ")";

	// Output node support as
	// internal node label
	// Do it reguardless of distances, because some consense
	// trees might have same tree strings
	// as other trees, we want to distinguish
	// between them in covSEARCH
	if (node->support) {
		char buf[16];
		sprintf(buf, "%d", node->support);
		ret = ret + buf;
	}

	ret = ret + "" + distance;	
	ret = ret;	
}


SCNode* 
SCTree::GetLeastSubtree(SCNode *node)
{
	if (node->IsLeaf())
		return node;

	SCNode* temp = node->children[0];
	bool allInts = true;
	unsigned int numChildren = node->NumChildren();

	for (unsigned int i = 0; i < numChildren; i++) {
		if (node->children[i]->IsLeaf()) {
			allInts = false;
			if (!temp->IsLeaf() || (temp->IsLeaf() && node->children[i]->name < temp->name))
				temp = node->children[i];
		}
	}

	if (allInts) {
		vector<SCNode*> kids;
		kids.resize(numChildren);

		for (unsigned int i = 0; i < numChildren; i++)
			kids[i] = GetLeastSubtree(node->children[i]);

		temp = kids[0];
		for (unsigned int i = 0 ;  i < numChildren; i++)
			if (kids[i]->name < temp->name)
				temp = kids[i];
	}

	return temp;
}

void 
SCTree::DeleteAllNodes()
{
	for (unsigned int i=0; i<nodelist.size(); ++i) {
		if (nodelist[i] != NULL) {
			delete nodelist[i];
			nodelist[i] = NULL;
		}
	}
}

// eof
