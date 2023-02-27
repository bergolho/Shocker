#include "node.h"

Node::Node (const uint32_t id)
{
    this->id = id;
    this->pos[0] = 0;
    this->pos[1] = 0;
    this->pos[2] = 0;
    this->is_active = false;
}

Node::Node (const uint32_t id, const double pos[])
{
    this->id = id;
    this->is_active = false;
    memcpy(this->pos,pos,sizeof(double)*3);
}

Node::Node (const uint32_t id, const double pos[], const bool is_active)
{
    this->id = id;
    this->is_active = is_active;
    memcpy(this->pos,pos,sizeof(double)*3);
}

Node::Node (Node *input)
{
    this->id = input->id;
    this->is_active = input->is_active;
    memcpy(this->pos,input->pos,sizeof(double)*3);
}

Node* Node::copy ()
{
    Node *result = new Node(); 
    result->id = this->id;
    result->is_active = this->is_active;
    memcpy(result->pos,this->pos,sizeof(double)*3); 
    return result;
}

Node::~Node ()
{
    
}

void Node::print ()
{ 
    printf("[point] id = %u || pos = (%g %g %g) || is_pmj = %d\n",this->id,this->pos[0],this->pos[1],this->pos[2],static_cast<int>(this->is_active));
}

void eliminate_node_from_list (std::vector<Node*> &n_list, Node *n)
{
    n_list.erase(n_list.begin() + n->id);
    order_node_list(n_list);
    delete n;
    n = nullptr;
}

void order_node_list (std::vector<Node*> &n_list)
{
    uint32_t counter = 0;
    for (uint32_t i = 0; i < n_list.size(); i++)
    {
        n_list[i]->id = counter;
        counter++;
    }
}

void print_node_list (std::vector<Node*> n_list)
{
    for (uint32_t i = 0; i < n_list.size(); i++)
        n_list[i]->print();
}

Node* search_node (std::vector<Node*> n_list, const uint32_t index)
{
    return n_list[index];
}

void copy_node_list (std::vector<Node*> in, std::vector<Node*> &out)
{
    for (uint32_t i = 0; i < in.size(); i++)
    {
        Node *n = new Node(in[i]);
        out.push_back(n);
    }
}

void concatenate_node_list (std::vector<Node*> in, std::vector<Node*> &out)
{
    uint32_t offset = out.size();
    for (uint32_t i = 0; i < in.size(); i++)
    {
        Node *n = new Node(in[i]);
        n->id += offset;
        out.push_back(n);
    }
}