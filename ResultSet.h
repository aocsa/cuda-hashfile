#pragma once
#ifndef __ResultSet_H
#define __ResultSet_H

#include <vector>
using std::vector;
 
 
template < class ObjectType >
class ResultSetPair{
   public:
      /**
      * Type of the object.
      */
      typedef ObjectType tObject;

      /**
      * Creates a new pair Object/Distance.
      *
      * <P>This instance will claim the ownership of the object. In other words,
      * it will dispose the object when it is no loger necessary.
      *
      * @param obj The result object.
      * @param distance The distance from the sample.
      */
      ResultSetPair(const tObject& obj, double distance){

         this->Object = obj;
         //this->JoinedObject = NULL;
         this->Distance = distance;
      }//end ResultSetPair

      /**
      * Creates a new pair Object/Distance with a JoinedObject.
      *
      * <P>This instance will claim the ownership of the object. In other words,
      * it will dispose the object when it is no loger necessary.
      *
      * @param obj The result object.
      * @param joined obj The result joined object.
      * @param distance The distance from the sample.
      */
      ResultSetPair(const tObject& obj, const tObject& joinedobj, double distance){

         this->Object = obj;
         this->JoinedObject = joinedobj;
         this->Distance = distance;
      }//end ResultSetPair

      /**
      * Disposes this instance and release all associated resources including
      * the result object.
      */
      ~ResultSetPair(){

      }//end ~ResultSetPair

      /**
      * This method returns the object.
      */
      const tObject GetObject(){
         return Object;
      }//end GetObject

      /**
      * This method returns the object.
      */
      const tObject GetJoinedObject(){
         return JoinedObject;
      }//end GetObject

      /**
      * This method gets the distance.
      */
      double GetDistance(){
         return Distance;
      }//end GetDistance
   public:
      /**
      * The object.
      */
      tObject Object;

      /**
      * The joined object.
      */
      tObject JoinedObject;

      /**
      * The distance from the query object.
      */
      double Distance;
};//end ResultSetPair

//----------------------------------------------------------------------------
// Class template ResultSet
//----------------------------------------------------------------------------
/**
* This class implements a query result. It will hold a set of pairs
* Object/Distance ordered by distance which constitutes the answer
* to a query. All query methods of all metric trees implemented by
* this library will return instances of this class.
*
* <P>In nearest neigbour queries and k-range queries, the result set may contain
* more than k results. It means that the greatest distance from sample to result
* has more than 1 object. In such cases, all objects whose distance from sample
* is equal to GetMaximumDistance() constitute the draw list.
*
* <P>As an extra capability, it can store information about the query but they
* will only be available if the query method supports this feature (this is an
* optional capability). See SetQueryInfo() for more details.
*
* <P> It also performs basic operations that allows the construction of
* the result set by the query procedures.
*
* @author Fabio Jun Takada Chino (chino@icmc.usp.br)
* @author Marcos Rodrigues Vieira (mrvieira@icmc.usp.br)
* @author Adriano Siqueira Arantes (arantes@icmc.usp.br)
* @version 1.1
* @ingroup struct
*/
template < class ObjectType >
class ResultSet{
   public:
      /**
      * Define constants to query type.
      */
      enum tQueryType{
         /**
         * Indicates that there is no information about the query.
         */
         UNKNOWN = 0,

         /**
         * This is a result of a range query.
         */
         RANGEQUERY = 1,

         /**
         * This is a result of a k nearest neighbour query.
         */
         KNEARESTQUERY = 2,

         /**
         * This is a result of a estimated k nearest neighbour query.
         */
         ESTIMATEKNEARESTQUERY = 3,

         /**
         * This is a result of a ring query.
         */
         RINGQUERY = 4,

         /**
         * This is a result of a k ring query.
         */
         KRINGQUERY = 5,

         /**
         * This is a result of a k AND range query.
         */
         KANDRANGEQUERY = 6,

         /**
         * This is a result of a k OR range query.
         */
         KORRANGEQUERY = 7,

         /**
         * This is a result of a Point query.
         */
         POINTQUERY = 8,

         /**
         * This is a result of a Reverse of Range query.
         */
         REVERSEDRANGEQUERY = 9,

         /**
         * This is a result of a K-Farthest query.
         */
         KFARTHESTQUERY = 10
      };//end tQueryType

      /**
      * Type of the object.
      */
      typedef ObjectType tObject;

      /**
      * This type defines the ResultSet Pair used by this class.
      */
      typedef ResultSetPair< ObjectType > tPair;

      /**
      * This method will create a new instance of this class. The parameter hint
      * is used to prepare this instance to hold at least <i>hint</i> results
      * (it is not a upper bound limit).
      *
      * @param hint The projected number of results (default = 1).
      */
      ResultSet(unsigned int hint = 1){
         // Reserve results
         Pairs.reserve(hint);
         // No info
         SetQueryInfo();
      }//end ResultSet

      /**
      * This method disposes this instance and releases all allocated resources.
      */
      virtual ~ResultSet();

      /**
      * This operator allow the access to a pair.
      */
      tPair & operator [] (unsigned int idx){
         return (*Pairs[idx]);
      }//end operator []

	  tPair & Get(unsigned int idx){
		  return (*Pairs[idx]);
	  }


      /**
      * This method returns the number of entries in this result.
      */
      unsigned int GetNumOfEntries(){
         return Pairs.size();
      }//end GetNumOfEntries

      /**
      * This method adds a pair Object/Distance to this result list.
      *
      * @param obj The object.
      * @param distance The distance from the sample object.
      * @warning There is no duplicate pair checking. All pairs will be added.
      */
      void AddPair(const tObject & obj, double distance){
         unsigned int pos;

         pos = this->Find(distance);
		 
         Pairs.insert(Pairs.begin() + pos, new tPair(obj, distance));
      }//end AddPair

      /**
      * This method adds a joined pair Object/JoinedObject/distance
      * to this result list.
      *
      * @param obj The object.
      * @param joinedObj The joined object.
      * @param distance The distance from the sample object.
      * @warning There is no duplicate pair checking. All pairs will be added.
      */
      void AddJoinedPair(const tObject& obj, const tObject& joinedObj, double distance){
         unsigned int pos;

         pos = this->Find(distance);
         Pairs.insert(Pairs.begin() + pos, new tPair(obj, joinedObj, distance));
      }//end AddJoinedPair


	  bool Contain(tObject* obj) {
		for (int i = 0; i < GetNumOfEntries(); i++) {
			if( Get(i).Object->IsEqual(obj)) {
				return true;
			}
		}
		return false;
	  }
      /**
      * This method will remove the last object from this result list.
      */
      void RemoveLast(){

         if (Pairs.size() > 0){
            if (Pairs[Pairs.size() - 1] != NULL){
               delete Pairs[Pairs.size() - 1];
            }//end if
            Pairs.pop_back();
         }//end if
      }//end RemoveLast

      /**
      * This method will remove the first object from this result list.
      */
      void RemoveFirst(){

         if (Pairs.size() > 0){
            if (Pairs[0] != NULL){
               delete Pairs[0];
            }//end if
            Pairs.erase(Pairs.begin());
         }//end if
      }//end RemoveFirst

      /**
      * This method returns the minimum distance of the objects in this result
      * list. If this result is empty, it will return a negative value.
      */
      double GetMinimumDistance(){

         if (Pairs.size() > 0){
            return Pairs[0]->GetDistance();
         }else{
            return -1;
         }//end if
      }//end GetMinimumDistance

      /**
      * This method returns the maximum distance of the objects in this result
      * list. If this result is empty, it will return a negative value.
      */
      double GetMaximumDistance(){

         if (Pairs.size() > 0){
            return Pairs[Pairs.size() - 1]->GetDistance();
         }else{
            return -1;
         }//end if
      }//end GetMaximumDistance

      /**
      * This method will cut out undesired objects. It may be used
      * by k-Nearest Neighbour queries to control the number of results.
      *
      * <P>This implementation also treat...
      *
      * @param limit The desired number of results.
      * @todo Review of this explanation.
      */
      void Cut(unsigned int limit);

      /**
      * This method will cut out undesired objects. It may be used
      * by k-Farthest Neighbour queries to control the number of results.
      *
      * <P>This implementation also treat...
      *
      * @param limit The desired number of results.
      * @todo Review of this explanation.
      */
      void CutFirst(unsigned int limit);

      /**
      * Adds information about the query. It is used by Query methods to add
      * information about the query. Since it is optional, not all results will
      * provide meaningful information about it.
      *
      * @param sample The sample object (a copy of it).
      * @param querytype The query type (UNKNOWN, RANGEQUERY, NEARESTQUERY,
      * KANDRANGEQUERY, KORRANGEQUERY, CROWNQUERY, KANDRANGEQUERYESTIMATE or
      * KORRANGEQUERYESTIMATE)
      * @param k The value of k (if it makes sence).
      * @param radius The value of radius (if it makes sence).
      * @param innerRadius The value of inner radius (if it makes sence).
      * @param tie The tie list. Default false;
      *
      * @warning It do not changes the behavior of this result.
      * @see GetQueryType()
      * @see GetK()
      * @see GetRadius()
      * @see GetSample()
      * @see GetTie()
      */
      void SetQueryInfo(const tObject& sample  = tObject(), int querytype = UNKNOWN,
                        int k = 0, double radius = 0.0,
                        double innerRadius = 0.0, bool tie = false){
         QueryType = querytype;
         K = k;
         Radius = radius;
         InnerRadius = innerRadius;
         Sample = sample;
         Tie = tie;
      }//end SetQueryInfo

      /**
      * Adds information about the query. It is used by Query methods to add
      * information about the query. Since it is optional, not all results will
      * provide meaningful information about it.
      *
      * @param sample The sample object (a copy of it).
      * @param querytype The query type (UNKNOWN, RANGEQUERY, KRANGEQUERY or
      * KNEARESTQUERY)
      * @param k The value of k (if it makes sence).
      * @param radius The value of radius (if it makes sence).
      * @param tie The tie list. Default false;
      *
      * @warning It do not changes the behavior of this result.
      * @see GetQueryType()
      * @see GetK()
      * @see GetRadius()
      * @see GetSample()
      * @see GetTie()
      */
      void SetQueryInfo(const tObject & sample , int querytype,
                        int k, double radius, bool tie){

			
			this->SetQueryInfo(sample, querytype, k, radius, 0.0, tie);
      }//end SetQueryInfo

      /**
      * Gets the information about query type. I may assume the values
      * UNKNOWN, RANGEQUERY, KANDRANGEQUERY, KORRANGEQUERY or KNEARESTQUERY.
      */
      int GetQueryType(){
         return QueryType;
      }//end GetQueryType

      /**
      * Gets information about k. It makes sense only for KANDRANGEQUERY, KORRANGEQUERY
      * and KNEARESTQUERY.
      */
      unsigned int GetK(){
         return K;
      }//end GetK

      /**
      * Gets information about radius. It makes sense only for RANGEQUERY and
      * KRANGEQUERY.
      */
      double GetRadius(){
         return Radius;
      }//end GetRadius

      /**
      * Gets information about inner radius. It makes sense only for CROWNQUERY.
      */
      double GetInnerRadius(){
         return InnerRadius;
      }//end GetRadius

      /**
      * Gets the sample object if it is available. Since it is an optional
      * information it may not be available.
      *
      * @return The sample object or NULL if it is not available.
      */
      const tObject& GetSample(){
         return Sample;
      }//end GetSample

      /**
      * Gets the tie list to the query. It is an optional
      * information.
      *
      * @return True if there is tie list or False otherwise.
      */
      bool GetTie(){
         return Tie;
      }//end GetTie

      /**
      * This method tests if two results are equal.
      *
      * @param r1 The second result to be test.
      * @return True if is equal, False otherwise.
      */
      bool IsEqual(ResultSet * r1);

      /**
      * This method tests the similarity between two results .
      *
      * @param r1 The second result to be test.
      * @return the percent-similarity of the two results.
      */
      double Precision(ResultSet * r1);


      /**
      * This method implements the intersection operator between two results.
      * @param r1 The first result set.
      * @param r2 The second result set.
      */
      void Intersection(ResultSet * result1, ResultSet * result2);

      /**
      * This method implements the union operator between two results.
      * @param r1 The first result set.
      * @param r2 The second result set.
      */
      void Union(ResultSet * result1, ResultSet * result2);

   private:
      /**
      * The vector of pairs.
      */
      vector < tPair * > Pairs;

      /**
      * Information about QueryType (Optional).
      */
      int QueryType;

      /**
      * Information about k for KANDRANGEQUERY, KORRANGEQUERY and
      * KNEARESTQUERY (Optional).
      */
      unsigned int K;

      /**
      * Information about tie list for KANDRANGEQUERY, KORRANGEQUERY
      * and KNEARESTQUERY .
      */
      bool Tie;

      /**
      * Information about radius for KANDRANGEQUERY, KORRANGEQUERY and
      * RANGEQUERY (Optional).
      */
      double Radius;

      /**
      * Information about inner radius for CROWNQUERY (Optional).
      */
      double InnerRadius;

      /**
      * Sample object (Optional).
      */
      tObject Sample;

      /**
      * This method locates the insertion position of an object.
      *
      * @param distance The desired distance.
      * @todo This code needs optimizations. I suggest a binary search
      * implementation.
      */
      unsigned int Find(double distance);
      
};//end ResultSet

/**********************************************************************
* GBDI Arboretum - Copyright (c) 2002-2004 GBDI-ICMC-USP
*
*                           Homepage: http://gbdi.icmc.usp.br/arboretum
**********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2004 Grupo de Bases de Dados e Imagens, Instituto de
 * Ciências Matemáticas e de Computação, University of São Paulo -
 * Brazil (the Databases and Image Group - Intitute of Matematical and 
 * Computer Sciences).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *       "This product includes software developed by Grupo de Bases
 *        de Dados e Imagens, Instituto de Ciências Matemáticas e de
 *        Computação, University of São Paulo - Brazil (the Databases 
 *        and Image Group - Intitute of Matematical and Computer 
 *        Sciences)"
 *
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names of the research group, institute, university, authors
 *    and collaborators must not be used to endorse or promote products
 *    derived from this software without prior written permission.
 *
 * 5. The names of products derived from this software may not contain
 *    the name of research group, institute or university, without prior
 *    written permission of the authors of this software.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OF THIS SOFTWARE OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * ====================================================================
 *                                            http://gbdi.icmc.usp.br/
 */
/**
* @file
*
* This file implements the classes ResultSet and ResultSets.
*
* @version 1.0
* $Revision: 1.15 $
* $Date: 2004/12/13 13:24:33 $
* $Author: marcos $
* @author Fabio Jun Takada Chino (chino@icmc.usp.br)
* @author Marcos Rodrigues Vieira (mrvieira@icmc.usp.br)
*/
// Copyright (c) 2002/2003 GBDI-ICMC-USP

//----------------------------------------------------------------------------
// Class template ResultSet
//----------------------------------------------------------------------------
template < class ObjectType >
ResultSet<ObjectType>::~ResultSet(){
   unsigned int i;

   for (i = 0; i < Pairs.size(); i++){
      if (Pairs[i] != NULL){
         //delete Pairs[i];//@todo: 
      }//end if
   }//end for

}//end ResultSet<ObjectType>::~ResultSet

//----------------------------------------------------------------------------
template < class ObjectType >
void ResultSet<ObjectType>::Cut(unsigned int limit){
   double max;
   bool stop;
   int i;

   // Will I do something ?
   if (GetNumOfEntries() > limit){
      if (Tie){ // if wants tie list
         // What is the max distance ?
         max = (*this)[limit - 1].GetDistance();

         // I'll cut out everybody which has distance greater than max.
         i = GetNumOfEntries() - 1;
         stop = i < limit;
         while (!stop){
            if ((*this)[i].GetDistance() > max){
               // Cut!
               RemoveLast();
               // The next to check is...
               i--;
               stop = (i < limit);
            }else{
               // Oops! I found someone who will not go out.
               stop = true;
            }//end if
         }//end while
      }else{
         RemoveLast();
      }//end if
   }//end if
}//end ResultSet<ObjectType>::Cut

//----------------------------------------------------------------------------
template < class ObjectType >
void ResultSet<ObjectType>::CutFirst(unsigned int limit){
   double min;
   bool stop;
   int idx;

   // Will I do something ?
   if (GetNumOfEntries() > limit){
      if (Tie){ // if wants tie list
         idx = GetNumOfEntries() - limit;
         // What is the min distance?
         min = (*this)[idx].GetDistance();
         // I'll cut out everybody which has distance lesser than min.
         stop = ((GetNumOfEntries() < limit) || (idx < 0));
         while (!stop){
            if ((*this)[idx].GetDistance() < min){
               // Cut!
               RemoveFirst();
            }//end if
            // The next to check is...
            idx--;
            stop = ((GetNumOfEntries() < limit) || (idx < 0));
         }//end while
      }else{
         RemoveFirst();
      }//end if
   }//end if
}//end ResultSet<ObjectType>::CutFirst

//----------------------------------------------------------------------------
template < class ObjectType >
unsigned int ResultSet<ObjectType>::Find(double distance){
   bool stop;
   unsigned int idx;

   idx = 0;
   stop = (idx >= Pairs.size());
   while (!stop){
      if (Pairs[idx]->GetDistance() < distance){
         idx++;
         stop = (idx >= Pairs.size());
      }else{
         stop = true;
      }//end if
   }//end while

   return idx;
}//end ResultSet<ObjectType>::Find

//----------------------------------------------------------------------------
template < class ObjectType >
bool ResultSet<ObjectType>::IsEqual(ResultSet * r1){
   ObjectType * tmp;
   bool result, result2;
   unsigned int i, j;
   int numObj1, numObj2;

   numObj1 = this->GetNumOfEntries();
   numObj2 = r1->GetNumOfEntries();

   // the default answer.
   result = false;
   // test if two results have the same number of entries and maximum
   // distance.
   if ((this->GetMaximumDistance() == r1->GetMaximumDistance()) && (numObj1 == numObj2)){
      // if there are, test one with the other.
      result = true;
      i = 0;
      // for each object in this class.
      while ((i < numObj1) && result) {
         // set the variables.
         result2 = false;
         // test starting with the first object. 
         j = 0;
         // for each object in the r1 set check until find a equal object.
         while ((j < numObj2) && (!result2)) {
            // check the distance first between the two objects.
            if (Pairs[i]->GetDistance() == r1->Pairs[j]->GetDistance()) {
               // the distance is equal, now test the object.
               tmp = r1->Pairs[j]->GetObject()->Clone();
               // set the result2 with the equality of the two objects.
               result2 = Pairs[i]->GetObject()->IsEqual(tmp);
               // delete the object's copy.
               delete tmp;
            }//end if
            // increment the couter of the second set.
            j++;
         }//end while
         // if the object in the first set was not in the second set, then
         // result will be false, otherwise true.
         result = result && result2;
         // increment the couter of the first set.
         i++;
      }//end while
   }//end if
   // return the result.
   return result;
}//end ResultSet<ObjectType>::IsEqual

//----------------------------------------------------------------------------
template < class ObjectType >
void ResultSet<ObjectType>::Intersection(ResultSet * result1, ResultSet * result2){
   bool result = false;
   ObjectType * tmpObj1, * tmpObj2;
   unsigned int idx, i;
   int numObj1, numObj2;
   double distance;

   numObj1 = result1->GetNumOfEntries();
   numObj2 = result2->GetNumOfEntries();

   if ((numObj1 != 0) && (numObj2 != 0)){
      for (idx = 0; idx < numObj1; idx++){

         // get the object from result1.
         tmpObj1 = (ObjectType *)result1->Pairs[idx]->GetObject();
         distance = result1->Pairs[idx]->GetDistance();
         i = 0;
         // check if the object from result1 is in result2.
         do{
            tmpObj2 = (ObjectType *)result2->Pairs[i]->GetObject();
            // is it equal to object1?
            result = tmpObj2->IsEqual(tmpObj1);
            // store if the two objects are equal.
            if (result){
              this->AddPair(tmpObj1->Clone(), distance);
            }//end if
            i++;
         }while ((i < numObj2) && (!result));
      }//end for
   }//end if

}//end ResultSet<ObjectType>::Intersection

//----------------------------------------------------------------------------
template < class ObjectType >
void ResultSet<ObjectType>::Union(ResultSet * result1, ResultSet * result2){
   bool result = false;
   ObjectType * tmpObj1, * tmpObj2;
   unsigned int idx, i;
   int numObj1, numObj2;
   double distance;

   numObj1 = result1->GetNumOfEntries();
   numObj2 = result2->GetNumOfEntries();

   if ((numObj1 != 0) && (numObj2 != 0)){
      // put all objects in result1 in unionResult.
      for (idx = 0; idx < numObj1; idx++){
         tmpObj1 = (ObjectType *)result1->Pairs[idx]->GetObject();
         distance = result1->Pairs[idx]->GetDistance();
         this->AddPair(tmpObj1->Clone(), distance);
      }//end for

      // now put all objects in result2 that are not in result1.
      for (idx = 0; idx < numObj2; idx++){
         // put all the objects in result1 and put in unionResult.
         tmpObj2 = (ObjectType *)result2->Pairs[idx]->GetObject();
         // it is storage the distance from result2's representative.
         distance = result2->Pairs[idx]->GetDistance();
         // check if the tmpObj2 in result2 is in result1.
         i = 0;
         do{
            tmpObj1 = (ObjectType *)result1->Pairs[i]->GetObject();
            distance = result1->Pairs[i]->GetDistance();
            // is it equal to object1?
            result = tmpObj1->IsEqual(tmpObj2);
            i++;
         }while ((i < numObj1) && (!result));
         // if the object2 is not in unionResult put it in unionResult.
         if (!result){
           this->AddPair(tmpObj2->Clone(), distance);
         }//end if
      }//end for
   }else if (numObj1 != 0){
      for (idx = 0; idx < numObj1; idx++){
         tmpObj1 = (ObjectType *)result1->Pairs[idx]->GetObject();
         distance = result1->Pairs[idx]->GetDistance();
         this->AddPair(tmpObj1->Clone(), distance);
      }//end for
   }else{
      for (idx = 0; idx < numObj2; idx++){
         tmpObj2 = (ObjectType *)result2->Pairs[idx]->GetObject();
         distance = result2->Pairs[idx]->GetDistance();
         this->AddPair(tmpObj2->Clone(), distance);
      }//end for
   }//end if
}//end ResultSet<ObjectType>::Union

//----------------------------------------------------------------------------
template < class ObjectType >
double ResultSet<ObjectType>::Precision(ResultSet * r1){
   double result = -1;
   ObjectType * tmpObj1, * tmpObj2;
   int idx, i;
   int numObj1, numObj2;
   int equal_count = 0;

   numObj1 = this->GetNumOfEntries();
   numObj2 = r1->GetNumOfEntries();

   // test if two results have the same number of entries
   if ((numObj1 != 0) && (numObj2 != 0)){
      for (idx = 0; idx < numObj2; idx++){
         // get the object from r1.
         tmpObj2 = (ObjectType *)r1->Pairs[idx]->GetObject();
         i = 0;
         // check if the object from r1 is in "this".
         do{
            tmpObj1 = (ObjectType *)this->Pairs[i]->GetObject();
            // is it equal to object1?
            result = tmpObj1->IsEqual(tmpObj2);
            // if the two objects are equal, increment the "equal_count".
            if (result){
               equal_count++;;
            }//end if
            i++;
         }while (i < numObj1); //end do while
      }//end for
   }//end if
   result = ((double)equal_count / numObj2);
   // return the result.
   return result;
}//end ResultSet<ObjectType>::Precision

#endif //__ResultSet_H
