/** 
* @mainpage Leipzig Array Package 
*
* @author Mario Fasold
*
* @section intro Introduction
*
* [Insert some introductionary text here].
*
* @section install Installation
*
* Following packes need to be installed to compile Larpack
*  @li libgsl0-dev (depends on libgsl-0) 
*  @li libboost-dev
*  @li [sometimes included in libboost-dev] libboost-program-options (-dev)
*
* @section CodingConventions Coding Conventions
* 
* @li Use coding conventions from Wrox: Professional C++
* @li Name Variables as accurate as possible. Qualifiers are placed at
*     at the end of a variable name (e.g.revenueTotal). Exception of
*     this rule are numerations (firstProbe, lastProbe, currentProbe).
*     Use Count and
*     Index instead of number (of), e.g. probeCount (total) and 
*     probeIndex. 
* @li Arrays are named in plural, if the new aggragation has no 
*     proper name (e.g. vector<Point> graph)
* @li Alwys write namespaces explicit in function definitions in headers
*     AND cpp files, even if that namespace is in a using clause. (e.g.
*     always int f(std::string& s)). Reason: This allows copy and paster
*     from header files (where using is not allowed) and doxygen needs it!  
* @li Simple getter and setter routines are defined in the CPP file as well,
*     to not plug up the header files with doxygen documentations.
*     Header Files may only include a brief decription of the function.
*     Exception are some templated classes, where all functions are defined
*     in the header file due to the "separate compiling model" (see
*     http://www.parashift.com/c++-faq-lite/templates.html#faq-35.12).
* @li Use doxygen comments in JavaDoc Style(@ instead of \). Classes and their
*     member variables are describted in the header file (.hpp). Member 
*     functions are described in the .cpp files, both brief and full description.
*     The first sentence is used as brief description.
*     Mathematical formulas in latex style are allowed.
* 
* @li [Functional Programming] Use a clean Strategy Pattern with a virtual
*     base class using the operator (). When speed matters, use boost::function
*     to pass a pointer to the function, thus avoiding dynamic binding.
*
*
* Subdirectories
* -# source/Math Mathematical functions.
* -# source/MicroarrayData Functions and classes handling the microarray data.
* -# source/Utility Utility functions e.g.String tools.
*
*
*
*
* [Test]
* <hr>
* @section notes release.notes
* release.notes
* <hr>
* @section requirements requirements
* @ verbinclude requirements
* <hr> 
*
*/
