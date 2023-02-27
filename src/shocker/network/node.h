//
// Created by bergolho on 18/08/21.
//

#ifndef NODE_H
#define NODE_H

#include <cstdio>
#include <cstdlib>
#include <cstdbool>
#include <cstdint>
#include <cstring>
#include <vector>

class Node
{
public:
    uint32_t id;
    double pos[3];
    bool is_active;
public:
    Node () { };
    Node (const uint32_t id);
    Node (const uint32_t id, const double pos[]);
    Node (const uint32_t id, const double pos[], const bool is_active);
    Node (Node *input);
    ~Node ();
    inline void setId (const uint32_t id) { this->id = id; }
    inline void setPosition (const double pos[]) { memcpy(this->pos,pos,sizeof(double)*3); }
    inline void setActive (const bool flag) { this->is_active = flag; }
    inline void setNode (Node *in) { this->id = in->id; memcpy(this->pos,in->pos,sizeof(double)*3); this->is_active = in->is_active; }
    Node* copy ();
    void print ();
};

Node* search_node (std::vector<Node*> n_list, const uint32_t index);
void eliminate_node_from_list (std::vector<Node*> &n_list, Node *n);
void order_node_list (std::vector<Node*> &n_list);
void print_node_list (std::vector<Node*> n_list);
void copy_node_list (std::vector<Node*> in, std::vector<Node*> &out);
void concatenate_node_list (std::vector<Node*> in, std::vector<Node*> &out);

#endif