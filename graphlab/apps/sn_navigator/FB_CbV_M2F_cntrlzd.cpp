/* This software is distributed under the MIT license (see LICENSE file in the
 * root directory)
 *
 * Copyright (c) 2013 Amira A. Soliman
 */

#include <graphlab.hpp>
#include <graphlab/graph/distributed_graph.hpp>
#include <graphlab/ui/metrics_server.hpp>
#include <graphlab/macros_def.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <algorithm>    // std::reverse
/*
*#################################### Data Structure Part #####################################
*/

size_t VIEW_SIZE = 50;
size_t PAYLOAD_SIZE = 10;


//////////////////////////////////////////////////////////
//gossip message exchanged between social graph nodes
//////////////////////////////////////////////////////////

struct single_message {
    graphlab::vertex_id_type sender_id;
    graphlab::vertex_id_type reciever_id;  
    std::multimap<std::string, std::string>  att_values;
    std::vector<std::string> att_group;
    std::vector<size_t> freq;

    single_message(): att_values() { }

   void save(graphlab::oarchive& oarc) const {
                oarc <<sender_id<<reciever_id;
                size_t num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                }
               num = att_group.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_group.begin();
                       iter != att_group.end(); ++iter){
                               oarc << *iter;
                }

               num = freq.size();
                oarc << num;
                for(std::vector<size_t>::const_iterator iter = freq.begin();
                       iter != freq.end(); ++iter){
                               oarc << *iter;
                }

   }
   void load(graphlab::iarchive& iarc) {
          iarc>>sender_id>>reciever_id;
          size_t num,f;
          std::string vk,v;
          iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }
         iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk;
                  att_group.push_back(vk);
               }
        iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq.push_back(f);
               }
   }
};

// message object exchanged using ivertex program
struct data_message: public graphlab::IS_POD_TYPE {
  std::vector<single_message> msgs;
  double value;

  data_message():  value(1.0), msgs() { }

  //
  // returns a priority to the reception of the message (used by the asynchronous engines)
  //
  double priority() const { return std::fabs(value); }

  //
  // allow message merging across the network
  // GraphLab always merges messages destined to the same vertex.
  //
  data_message& operator+=(const data_message& other) {

          for(size_t i = 0;i < other.msgs.size(); ++i){
               msgs.push_back(other.msgs[i]);
          }

          return *this;
  }
  void save(graphlab::oarchive& oarc) const {

        size_t num_msgs = msgs.size();
        oarc << num_msgs;
        for(size_t k = 0;k < num_msgs; ++k){
                 msgs[k].save(oarc);
        }
  }

  void load(graphlab::iarchive& iarc) {

          size_t num_msgs = 0;
          msgs.clear();

          iarc >> num_msgs;

          for(size_t k = 0;k < num_msgs; ++k){
                  single_message msg;
                  msg.load(iarc);
                  msgs.push_back(msg);
          }
  }
};
//////////////////////////////////////////////////////////
//edge
struct edata {
  graphlab::edge_id_type weight;

  edata():weight(0) {
  }

  explicit edata(graphlab::edge_id_type w): weight(w) {
  }

  void save(graphlab::oarchive& oarc) const {
    oarc << weight;
  }

  void load(graphlab::iarchive& iarc) {
    iarc >> weight;
  }

};
////////I//////////////////////////////////////////////////
//vertex
struct vdata {

        // list of direct friends
                
        std::vector<graphlab::vertex_id_type> friends;
        std::vector<std::string> att_names;
        std::multimap<std::string, std::string>  att_values;
        std::multimap<std::string, std::vector<std::string> > att_v_n_freq;
        std::vector<std::string> att_group;
        std::vector<size_t> freq;
        std::vector<double> CAG_supp;
        std::vector<std::string> att_group_clstr;
        std::vector<size_t> freq_clstr;

     vdata() :
                friends(), att_names(), att_values(), att_v_n_freq(), att_group(),freq(), att_group_clstr(),freq_clstr(),CAG_supp() {}

      explicit vdata( std::multimap<std::string, std::string> a_v) :
               friends(),att_names(), att_values(a_v), att_v_n_freq(), att_group(),freq(), att_group_clstr(),freq_clstr() ,CAG_supp() {}

        void save(graphlab::oarchive& oarc) const {
               
                //store list of friends
                
                size_t num = friends.size();
                oarc << num;
                for(size_t i = 0;i < num; ++i){
                        oarc << friends[i];
                }

                //CbIV
                num = att_names.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_names.begin();
                       iter != att_names.end(); ++iter){
                               oarc << *iter;
                }
                
                num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
                }
                num = att_group.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_group.begin();
                       iter != att_group.end(); ++iter){
                               oarc << *iter;
                }
                num = freq.size();
                oarc << num;
                for(std::vector<size_t>::const_iterator iter = freq.begin();
                       iter != freq.end(); ++iter){
                               oarc << *iter;
                }
                num = CAG_supp.size();
                oarc << num;
                for(std::vector<double>::const_iterator iter = CAG_supp.begin();
                       iter != CAG_supp.end(); ++iter){
                               oarc << *iter;
                }

                num = att_group_clstr.size();
                oarc << num;
                for(std::vector<std::string>::const_iterator iter = att_group_clstr.begin();
                       iter != att_group_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = freq_clstr.size();
                oarc << num;
                for(std::vector<size_t>::const_iterator iter = freq_clstr.begin();
                       iter != freq_clstr.end(); ++iter){
                               oarc << *iter;
                }
                num = att_v_n_freq.size();
                oarc << num;
                for(std::multimap<std::string, std::vector<std::string> >::const_iterator iter = att_v_n_freq.begin();
                       iter != att_v_n_freq.end(); ++iter){
                               oarc << iter->first;
                               std::vector<std::string> vlist =  iter->second;
                               for(size_t vi =0; vi<vlist.size(); vi++)
                                   oarc << vlist[vi];;
                }
                



        }
     void load(graphlab::iarchive& iarc) {
               
               std::string v,vk;
               size_t num = 0;

                iarc >> num;
                for(size_t i = 0; i < num; i++) {
                  graphlab::vertex_id_type f_id = 0;
                  iarc >> f_id;
                  friends.push_back(f_id);
                }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk ;
                  att_names.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
                  //if(vk != "interest")
                  //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
               }

               iarc >> num;
               size_t f=0;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk;
                  att_group.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq.push_back(f);
               }
               double d;
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> d;
                  CAG_supp.push_back(d);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk;
                  att_group_clstr.push_back(vk);
               }

               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> f;
                  freq_clstr.push_back(f);
               }
               iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  std::vector<std::string> vs;
                  iarc >> vk ;
                  //Att1value
                  iarc >> v;
                  vs.push_back(v);
                  //Att2value
                  iarc >> v;
                  vs.push_back(v);
                  //Att3value
                  iarc >> v;
                  vs.push_back(v);

                  att_v_n_freq.insert(std::make_pair(vk,vs));
               }
        }

 };


///////////////////////////////global functions/////////////////////////////////

int find(std::vector<std::string> vec, std::string key){
    int index =-1;

    for(size_t i=0; i<vec.size(); i++){
       if(vec[i] == key){
          index =i; 
          break;
       }
    }
    return index;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//distributed graph definition
typedef graphlab::distributed_graph<vdata, edata> graph_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//#################################### Send data Items part #####################################
//

struct data_gatherer_directN {

        std::multimap<std::string, std::string> att_values;

        data_gatherer_directN () : att_values() {}
        void save(graphlab::oarchive& oarc) const {
            size_t num = att_values.size();
                oarc << num;
                for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin();
                       iter != att_values.end(); ++iter){
                               oarc << iter->first;
                               oarc << iter->second;
                      std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
                }

        }
        void load(graphlab::iarchive& iarc) {
             size_t num;
             std::string vk, v;
             att_values.clear();
             iarc >> num;
               for(size_t i = 0; i < num; i++) {
                  iarc >> vk >> v;
                  att_values.insert(std::make_pair(vk,v));
               }

        }
       data_gatherer_directN& operator+=(const data_gatherer_directN& other) {
              for(std::multimap<std::string, std::string>::const_iterator iter = other.att_values.begin();
                  iter != other.att_values.end(); ++iter){
                  att_values.insert(std::make_pair(iter->first,iter->second));
                   //if(iter->first!="interest")
                   //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }
              return *this;
         }
};

class collect_attributes:public graphlab::ivertex_program<graph_type,data_gatherer_directN,data_message>, public graphlab::IS_POD_TYPE {

public:
 
       std::multimap<std::string, std::string> att_values;
       std::vector<single_message> msg2frwd;
       std::vector<std::string> att_group_clstr;
       std::vector<size_t> freq_clstr;

        edge_dir_type gather_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::ALL_EDGES;
        }

        edge_dir_type  scatter_edges(icontext_type& context,
                              const vertex_type& vertex) const {
                return graphlab::NO_EDGES;
        }

        data_gatherer_directN gather(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
          data_gatherer_directN att_n_v;
          //std::cout << vertex.id() << ", " << vertex.data().friends.size() <<"\n";
          
          
             /*for(std::multimap<std::string, std::string>::const_iterator iter = vertex.data().att_values.begin();
                       iter != vertex.data().att_values.end(); ++iter) {
                       att_n_v.att_values.insert(std::make_pair(iter->first,iter->second));
                       //if( iter->first =="education")
                       //std::cout <<"key: "<<iter->first<<", val: "<<iter->second<< "\n";
              }
              
             if(edge.target().id() == vertex.id() && (edge.source().data().att_values.size() >0))               
               for(std::multimap<std::string, std::string>::iterator iter = edge.source().data().att_values.begin(); 
                   iter!=edge.source().data().att_values.end(); iter++){
                   std::string vk = iter->first;
                   std::string v = iter->second;
                   att_n_v.att_values.insert(std::make_pair(vk, v));
                   //if( vk !="interest")
                   //std::cout <<"key: "<<vk<<", val: "<<v<< "\n";
              }
         
           */
        
          return att_n_v;
        }

        void scatter(icontext_type& context, const vertex_type& vertex,
               edge_type& edge) const {
        }


   //====================================================================================

   void init(icontext_type& context,
                const vertex_type& vertex, const data_message &msg) {

          if(context.iteration() == 1 ){
          for(size_t i=0; i<msg.msgs.size() ; i++){
             if(vertex.id() == msg.msgs[i].reciever_id){
                 single_message reply =  msg.msgs[i];
                 for(std::multimap<std::string, std::string>::iterator iter = reply.att_values.begin();
                   iter!=reply.att_values.end(); iter++){
                   std::string vk = iter->first;
                   std::string v = iter->second;
                   att_values.insert(std::make_pair(vk, v));
             }
          }
        }
      } 
   }

   void apply(icontext_type& context, vertex_type& vertex,const gather_type& total) {
         
    
    if(context.iteration()==0){
       single_message reply;
       reply.sender_id = vertex.id();
       reply.reciever_id = 595587084;
       reply.att_values = vertex.data().att_values;
       data_message msg;
       msg.msgs.push_back(reply);
      context.signal_vid(595587084, msg);
    }
    else{
       std::vector<std::string> att_names;
       att_names.push_back("gndr");
       att_names.push_back("fname");
       att_names.push_back("education");
       att_names.push_back("job");
       att_names.push_back("employer");
       att_names.push_back("interest");
       att_names.push_back("home_cntry");
       att_names.push_back("currloc_cntry");
      
           for(std::multimap<std::string, std::string>::const_iterator iter = att_values.begin(); iter!= att_values.end(); ++iter)
               vertex.data().att_values.insert(std::make_pair(iter->first, iter->second));
           //calculting the frequencies using node's att_values
           for(size_t i=0;i<att_names.size()-1; i++){
               std::string a_n1 = att_names[i];
               std::vector<std::string> triple;
               for(size_t j=i+1; j<att_names.size(); j++){
                  std::string a_n2 = att_names[j];   
                  if(a_n1 != a_n2){
                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator > 
                       iter_n1 = vertex.data().att_values.equal_range(a_n1);
                  std::pair<std::multimap<std::string, std::string>::iterator,std::multimap<std::string, std::string>::iterator >
                       iter_n2 = vertex.data().att_values.equal_range(a_n2);
                  std::vector<std::string> val_n1, val_n2;
                  for (std::multimap<std::string, std::string>::iterator it=iter_n1.first; it!=iter_n1.second; ++it)
                      val_n1.push_back(it->second);
                  for (std::multimap<std::string, std::string>::iterator it=iter_n2.first; it!=iter_n2.second; ++it)
                      val_n2.push_back(it->second);
                  std::string value1, value2;
                  size_t tot_tran_size =  val_n1.size() + val_n2.size();
                  size_t freq = 1;
                  size_t len = val_n1.size();
                  if(len > val_n2.size())
                     len = val_n2.size();
                  for(size_t vali =0; vali<len; vali ++){
                     value1 = val_n1[vali];
                     value2 = val_n2[vali];
                     //start to find repetition and increase counter//check if counting is done before for same values
                     std::multimap<std::string,std::vector<std::string> >::iterator fit = vertex.data().att_v_n_freq.find(value1+ "&" + value2);
                     if(fit == vertex.data().att_v_n_freq.end()){
                        triple.clear();
                        triple.push_back(a_n1);
                        triple.push_back(a_n2);     
                        for(size_t index = vali +1; index<len; index++)
                           if(value1 == val_n1[index] && value2 == val_n2[index])
                              freq ++; 
                        std::stringstream convert;
                        convert << freq;
                        triple.push_back(convert.str());
                        convert << tot_tran_size;
                        triple.push_back(convert.str());
                        vertex.data().att_v_n_freq.insert(std::make_pair(value1+ "&" + value2, triple));
                     }  
                  }
                  }
               }
            }
           //add final aggregation on discovered CAG
            for(std::multimap<std::string, std::vector<std::string> >::iterator iter = vertex.data().att_v_n_freq.begin();
                       iter != vertex.data().att_v_n_freq.end(); ++iter){
                               std::vector<std::string> vlist =  iter->second;
                               std::string key = vlist[0] + "&" +  vlist[1];
                               size_t freq = atoi(vlist[2].c_str());
                               double supp =  double(atoi(vlist[2].c_str()))/double(atoi(vlist[3].c_str()));
                               
                               int index = find(vertex.data().att_group, key);
                               if(index ==-1){
                                  vertex.data().att_group.push_back(key);
                                  vertex.data().freq.push_back(freq);
                                  vertex.data().CAG_supp.push_back(supp);
                               }
                               else{
                                   vertex.data().freq[index] += freq;
                                   vertex.data().CAG_supp[index] += supp;
                               }
              
                }
           }   
   }
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//############################ line parser to load data from file ##################################################//


bool line_parser(graph_type& graph,
                   const std::string& filename,
                   const std::string& textline)
 {
    if(textline.empty()) return true;
 
    std::stringstream strm(textline);

    graphlab::vertex_id_type vid;
    std::vector<graphlab::vertex_id_type> frnds;
    std::string attribute, token, token_sub, token_list;
    bool  loadmiss = false;
    std::vector<std::string> att_names;
    std::multimap<std::string, std::string>  att_values;

/////////////////////////////////////////////////Load data    
   strm >> vid;

    strm >> attribute;

    //graph.add_vertex(vid);

    if(attribute != "NULL") {
       
       std::istringstream ss(attribute);
       while(std::getline(ss, token, ',')) {
          size_t fID = atoi(token.c_str());
          if(fID != vid){
             frnds.push_back(fID);
             //graph.add_edge(fID, vid);
          }
       }
    }


   strm >> attribute;
   //std::cout << "fname " << attribute << "\n";
   if(attribute != "NULL")
   att_values.insert(std::make_pair("fname",attribute));

   strm >> attribute;
   //std::cout << "gndr " << attribute << "\n";
   if(attribute != "NULL")
   att_values.insert(std::make_pair("gndr",attribute));

   //strm >> attribute;
   //att_values.insert(std::make_pair("home_city",attribute));

   strm >> attribute;
   //std::cout << "hc " << attribute << "\n";
   if(attribute != "NULL")
   att_values.insert(std::make_pair("home_cntry",attribute));

   //strm >> attribute;
   //att_values.insert(std::make_pair("currloc_city",attribute));

   strm >> attribute;
   //std::cout << "curc " << attribute << "\n";
   if(attribute != "NULL")
   att_values.insert(std::make_pair("currloc_cntry",attribute));

   strm >> attribute;
   if(attribute != "NULL") {
    std::istringstream ss1(attribute);
   
   //if(wht2load<=10 || more2load<=10){
 
   while(std::getline(ss1, token, ',')){
      std::istringstream ss12(token);
       while(std::getline(ss12, token_sub, '_'))
          att_values.insert(std::make_pair("job",token_sub));
    }
   }
   strm >> attribute;
   if(attribute != "NULL") {
   std::istringstream ss2(attribute);
   //if(wht2load<=10 || more2load<=10)
   while(std::getline(ss2, token, ',')){
       std::istringstream ss12(token);
       while(std::getline(ss12, token_sub, '_'))
	  att_values.insert(std::make_pair("employer",token_sub));
   }
   }
   //strm >> attribute;
   //std::istringstream ss9(attribute);
   //while(std::getline(ss9, token, ','))
      //att_values.insert(std::make_pair("job_city",token));

   //strm >> attribute;
   //std::istringstream ss3(attribute);
   //while(std::getline(ss3, token, ','))
      //att_values.insert(std::make_pair("job_cntry",token));

   strm >> attribute;
   if(attribute != "NULL") {
   std::istringstream ss4(attribute);
   
   //if((wht2load>10 && wht2load<=20) || (more2load>10 && more2load<=20)){
   
   while(std::getline(ss4, token, ',')){
      std::istringstream ss12(token);
       while(std::getline(ss12, token_sub, '_'))
      att_values.insert(std::make_pair("education",token_sub));
      /*int index =0;
      std::string cat;
      std::istringstream ss12(token);
       while(std::getline(ss12, token_sub, '$')){
         if(index ==0){
            cat = token_sub;
            index ++;
         }
         else{
          att_values.insert(std::make_pair("education",cat + "_" + token_sub));
          //std::cout << "edu " << cat + "_" + token_sub << "\n"; 
         }
       }*/
   }
  //}
  }
   //strm >> attribute;
   //std::istringstream ss5(attribute);
   //while(std::getline(ss5, token, ','))
      //att_values.insert(std::make_pair("edu_org",token));

   //strm >> attribute;
   //std::istringstream ss6(attribute);
   //while(std::getline(ss6, token, ','))
       //att_values.insert(std::make_pair("edu_city",token));

   //strm >> attribute;
   //std::istringstream ss7(attribute);
   //while(std::getline(ss7, token, ','))
      //att_values.insert(std::make_pair("edu_cntry",token));
   
   strm >> attribute;
   if(attribute!="NULL"){
   std::istringstream ss8(attribute);
   
   while(std::getline(ss8, token, ',')){
      std::istringstream ss12(token);
       int index=0;
       std::string cat="";
      while(std::getline(ss12, token_sub, '@')){
         if(index==0){
            cat= token_sub;
            index++;
         }
         else {
           std::istringstream ss(token_sub);
           while(std::getline(ss, token_list, '_'))
              att_values.insert(std::make_pair("interest",cat + "#" + token_list));
         }
      }
   }
   }
   
  //for(std::multimap<std::string, std::string>::iterator iter = att_values.begin(); iter!=att_values.end(); iter++)
           //std::cout << "key: " << iter->first << ", val: " << iter->second << "\n"; 


///////////////////////////////////////////Add node to the graph
    graph.add_vertex(vid, vdata(att_values));
    if(frnds.size()>0)
       //std::cout << vid << " , " << frnds.size() << "\n";
    for(int  i=0; i< frnds.size(); i++){
       graph.add_edge(frnds[i], vid);
       //graph.add_edge( vid, frnds[i]);
    }


   return true;
 }


//#################################  Graph Writer ######################################//
//#################################  Graph Writer ######################################//
struct graph_writer {
  std::string save_vertex(graph_type::vertex_type v) {


        std::stringstream strm;
        std::vector<std::string> att_vec;
        std::vector<size_t> att_freq;
         //strm << "============================================== \n";
         //strm <<  v.id() << "\t num of friends:" << v.data().friends.size() << "\n";
        if(v.id()==595587084){
         strm << "============================================== \n";
         strm <<  v.id() << "\t num of friends:" << v.data().friends.size() << "\n";
         for(size_t i=0; i< v.data().att_group.size(); i++)
             strm <<  v.data().att_group[i] <<"," <<  v.data().freq[i] <<","  << v.data().CAG_supp[i] << "\n";
         }
         /*strm << "************************ \n";
         for(std::multimap<std::string, std::vector<std::string> >::iterator iter = v.data().att_v_n_freq.begin();
                       iter != v.data().att_v_n_freq.end(); ++iter){
                   std::vector<std::string> vlist =  iter->second;
                   size_t freqv = atoi(vlist[2].c_str());
                   strm <<  iter->first  <<"," <<  freqv << "\n";
         }*/
         //strm << "************************ \n";
         
         //for(size_t i=0; i< v.data().att_group_clstr.size(); i++)
             //strm <<  v.data().att_group_clstr[i] <<"," <<  v.data().freq_clstr[i] << "\n";
        //}
    return strm.str();
  }
  std::string save_edge(graph_type::edge_type e) { return ""; }
};

//######################################## main function #########################################//
//######################################## main function #########################################//
  int main(int argc, char** argv) {

        graphlab::command_line_options clopts("Dosn: fill neighbors map");
        std::string file_name;
        clopts.attach_option("graph", file_name, "Graph input. reads all graphs matching prefix*");

        if(!clopts.parse(argc, argv)) return EXIT_FAILURE;
        if (file_name == "") {
                std::cout << "--graph is not optional\n";
                clopts.print_description();
                return EXIT_FAILURE;
        }

        // Initialize control plane using mpi
        graphlab::mpi_tools::init(argc, argv);
        graphlab::distributed_control dc;
        graphlab::launch_metric_server();

        // load graph
        graph_type graph(dc, clopts);      //graph_type graph(dc);
         
        graph.load(file_name, line_parser); //load_from_hdfs("hdfs://graphLabLarge1:154310/data/graph.txt", line_parser);
        //graph.load_format(file_name, "tsv");
        graph.finalize();
        dc.cout() << "Number of vertices: " << graph.num_vertices() << std::endl
              << "Number of edges:    " << graph.num_edges() << std::endl;

   
   size_t MAX_ITERATIONS=2;
   clopts.get_engine_args().set_option("max_iterations", MAX_ITERATIONS);
   graphlab::omni_engine<collect_attributes> engine2(dc, graph, "sync", clopts);
   engine2.signal_all();
   engine2.start();

   graph.save("FB_Data.txt", graph_writer(),
             false,    // do not gzip
             true,     // save vertices
             false);   // do not save edges
   
   graphlab::mpi_tools::finalize();



   return EXIT_SUCCESS;
  }



