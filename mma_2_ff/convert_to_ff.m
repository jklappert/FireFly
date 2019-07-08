(*//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================*)
convert[funs_, vars_, nthr_] :=
    Module[{stream,stream2,rem,tmplist,tlist,clist},	   
	   CreateFile["ff_conv/funs.hpp"];
	   CreateFile["Makefile"];

	   CreateFile["ff_conv/exec.cpp"];

	   stream2 = OpenWrite["ff_conv/exec.cpp"];
	   WriteString[stream2, "#include \"Reconstructor.hpp\"\n"];
	   WriteString[stream2, "#include \"funs.hpp\"\n\n"];
	   WriteString[stream2, "using namespace firefly;\n\n"];
           WriteString[stream2, "class BlackBoxUser : public BlackBoxBase{\npublic:\n  BlackBoxUser(){};\n"];
	   WriteString[stream2, "  virtual std::vector<FFInt> operator()(const std::vector<FFInt>& values){\n  }\n"];
	   WriteString[stream2, "  virtual void prime_changed(){\n  }\n};\n"];
	   WriteString[stream2, "int main() {\n"];
	   WriteString[stream2, "  BlackBoxUser bb;\n"];
	   WriteString[stream2, "  Reconstructor reconst("<>ToString[Length[vars]]<>","<>ToString[nthr]<>",bb);\n"];
	   WriteString[stream2, "  reconst.enable_scan();\n"];
	   WriteString[stream2, "  reconst.reconstruct();\n"];
	   WriteString[stream2, "  return 0;\n}\n\n"];
	   Close[stream2];
	   
	   stream = OpenWrite["ff_conv/funs.hpp"];
	   WriteString[stream, "#pragma once\n"];
	   WriteString[stream, "#include \"FFInt.hpp\"\n"];
	   WriteString[stream, "namespace firefly{\n"];

	   stream2 = OpenWrite["Makefile"];
	   WriteString[stream2, "CXX      := -c++\n"];
           WriteString[stream2, "CXXFLAGS := -O2 -std=c++11\n"];
	   WriteString[stream2, "LIBPATHS := -L$(FF_LIB_PATH) -L$(FLINT_LIB_PATH)\n"];
	   WriteString[stream2, "LDFLAGS  := -lfirefly -lflint -lgmp\n"];
	   WriteString[stream2, "BUILD    := ./build\n"];
	   WriteString[stream2, "OBJ_DIR  := $(BUILD)/objects\n"];
	   WriteString[stream2, "EXEC_DIR := $(BUILD)\n"];
	   WriteString[stream2, "TARGET   := exec\n"];
	   WriteString[stream2, "INCLUDE  := -Iff_conv/\n"];
	   WriteString[stream2, "SRC      := $(wildcard ff_conv/*.cpp)\n\n"];
	   WriteString[stream2, "OBJECTS  := $(patsubst ff_conv/%.cpp,$(OBJ_DIR)/%.o,$(SRC))\n\n"];
	   WriteString[stream2, "all: build $(EXEC_DIR)/$(TARGET)\n\n"];
	   WriteString[stream2, "$(OBJ_DIR)/%.o: ff_conv/%.cpp\n"];
	   WriteString[stream2, "\t@mkdir -p $(@D)\n"];
	   WriteString[stream2, "\t$(CXX) $(CXXFLAGS) -I$(FF_INC_DIR) $(INCLUDE) -o $@ -c $<\n\n"];
	   WriteString[stream2, "$(EXEC_DIR)/$(TARGET): $(OBJECTS)\n"];
	   WriteString[stream2, "\t@mkdir -p $(@D)\n"];
	   WriteString[stream2, "\t$(CXX) $(CXXFLAGS) -I$(FF_INC_DIR) $(INCLUDE) $(LIBPATHS) -o $(EXEC_DIR)/$(TARGET) $(OBJECTS) $(LDFLAGS)\n\n"];

	   WriteString[stream2, ".PHONY: all build clean\n\n"];
	   WriteString[stream2, "build:\n"];
	   WriteString[stream2, "\t@mkdir -p $(EXEC_DIR)\n"];
	   WriteString[stream2, "\t@mkdir -p $(OBJ_DIR)\n\n"];
	   WriteString[stream2, "clean:\n"];
	   WriteString[stream2, "\t-@rm -rvf $(OBJ_DIR)/*\n"];
	   WriteString[stream2, "\t-@rm -rvf $(EXEC_DIR)/*\n\n"];
	   Close[stream2];

	   Do[
	       Print["Converting function "<>ToString[ii]<>"..."];
	       If[SameQ[Head[funs[[ii]]], String],
		  num = Expand[Numerator[ToExpression[funs[[ii]]]]];
		  den = Expand[Denominator[ToExpression[funs[[ii]]]]];
		  ,
		  num = Expand[Numerator[funs[[ii]]]];
                  den = Expand[Denominator[funs[[ii]]]];
		 ];

	       If[SameQ[Head[num], Plus],
		  tmplist = List@@Expand[num];
		  rem = QuotientRemainder[Length[num],500];
		  tlist = {rem[[2]]};
		  (* Get the splitting of the coefficients *)
		  kk = 1;
		  Do[
		      tlist = Join[tlist, {500}];,
		      {kk,1,rem[[1]]}
		    ];
		  (* Split the coefficients *)
		  clist = TakeList[tmplist, tlist];,
		  clist = {{num}};
		 ];

	       (* write terms to files *)
	       writeFuns[clist,vars,stream,"num",ii];

	       (* same procedure for denominator *)
	       If[SameQ[Head[den], Plus],
		  len = Length[den];
		  tmplist = List@@Expand[den];
		  rem = QuotientRemainder[len,500];
		  tlist = {rem[[2]]};

		  (* Get the splitting of the coefficients *)
		  Do[
		      tlist = Join[tlist, {500}];,
		      {kk,1,rem[[1]]}
		    ];
		  
		  (* Split the coefficients *)
		  clist = TakeList[tmplist, tlist];,
		  clist = {{den}};
		 ];

	       (* write terms to files *)
	       writeFuns[clist,vars,stream,"den",ii];,
	       {ii,1,Length[funs]}
	     ];

	   (* Write functions to header *)
	   Do[
               WriteString[stream, "  inline FFInt fun"<>ToString[ii]<>"(const std::vector<FFInt>& yis){return fun"<>ToString[ii]<>"num(yis) / fun"<>ToString[ii]<>"den(yis);};\n"];,
               {ii,1,Length[funs]}
             ];

	   WriteString[stream, "}\n"];
           Close[stream];
	  ];

writeFuns[coefs_,vars_,strm_,suff_,ii_] :=
    Module[{sum},
	   Do[
	       If[Length[coefs[[kk]]] < 200,
			 tmpterms = Simplify[Apply[Plus,coefs[[kk]]]];,
			 tmpterms = Collect[Apply[Plus,coefs[[kk]]],ToExpression[vars], Simplify];
			];
	       tmpterms = tmpterms /. Rational[a_,b_] :> FFInt["mpz_class(",ToString[a],")"]/FFInt["mpz_class(",ToString[b],")"] /; Length[IntegerDigits[a]] > 18 && Length[IntegerDigits[b]] > 18;
	       tmpterms = tmpterms /. Rational[a_,b_] :> FFInt["mpz_class(",ToString[a],")"]/FFInt[b] /; Length[IntegerDigits[a]] > 18;
	       tmpterms = tmpterms /. Rational[a_,b_] :> FFInt[a]/FFInt["mpz_class(",ToString[b],")"] /; Length[IntegerDigits[b]] > 18;
	       tmpterms = tmpterms /. Rational[a_,b_] :> FFInt[a]/FFInt[b];
	       tmpterms = tmpterms /. a_ :> FFInt["mpz_class(",ToString[a],")"] /; Length[IntegerDigits[a]] > 18(*IntegerQ[a]*);
	       (*	       tmpterms = tmpterms //. Power[a__,b_] :> pow[a,b] /; MemberQ[ToExpression[vars],a];*)
	       tmpterms = tmpterms //. Power[a__,b_] :> pow[a,b] /; b > 0;
(*	       tmpterms = tmpterms /. a_ :> FFInt[a] /; IntegerQ[a];*)
	       CreateFile["ff_conv/fun"<>ToString[ii]<>suff<>ToString[kk]<>".cpp"];
	       stream3 = OpenWrite["ff_conv/fun"<>ToString[ii]<>suff<>ToString[kk]<>".cpp"];
	       WriteString[stream3, "#include \"funs.hpp\"\n\n"];
               WriteString[stream3, "namespace firefly{\n\n"];
	       WriteString[stream3, "  FFInt fun"<>ToString[ii]<>suff<>ToString[kk]<>"(const std::vector<FFInt>& yis){\n"];
	       Do[
		   WriteString[stream3, "    FFInt "<>ToString[vars[[jj]]]<>" = yis["<>ToString[jj -1]<>"];\n"];,
                   {jj,1,Length[vars]}
                 ];

	       WriteString[stream3, "\n"];
               WriteString[stream3, "    FFInt "<>suff<>" = "<>ToString[tmpterms,CForm]<>";\n"];
               WriteString[stream3, "    return "<>suff<>";\n"];
               WriteString[stream3, "  }\n}"];
               Close[stream3];,
	       {kk,1,Length[coefs]}
	     ];

	   sum = "0";
	   Do[
	       WriteString[strm, "  FFInt fun"<>ToString[ii]<>suff<>ToString[kk]<>"(const std::vector<FFInt>& yis);\n"];
	       sum = sum<>" + fun"<>ToString[ii]<>suff<>ToString[kk]<>"(yis)";
	       ,
	       {kk,1,Length[coefs]}
	     ];
	   WriteString[strm, "  inline FFInt fun"<>ToString[ii]<>suff<>"(const std::vector<FFInt>& yis){return "<>sum<>";};\n"];
	  ];

If[StringLength[functions] > 0,
   Print["Converting "<>functions<>" with variables "<>variables<>" and "<>ToString[nthreads]<>" threads."];
   funs = ToExpression[ReadString[functions]];
   vars = ToExpression[ReadString[variables]];
   convert[funs, vars, nthreads];
  ];

