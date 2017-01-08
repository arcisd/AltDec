
#--------------------------------------------------------------------------------------------------------------------------------------------

# By Diego ARCIS / 2013. / to compute heads, tales, and right alternating decompositions of elements in Artin groups.

# Use GAP3 with Chevie Package 
# GAP3 (including Chevie Package) may be obtained from "https://webusers.imj-prg.fr/~jean.michel/gap3/".

# B:=Braid(W); / braid group (Artin group) of type "W", where "W" is one of the following Coxeter groups:

# W:=CoxeterGroup("A",k);   # (k>0 generators)
# W:=CoxeterGroup("B",k);   # (k>1 generators)
# W:=CoxeterGroup("D",k);   # (k>3 generators)
# W:=CoxeterGroup("I",2,t); # (t>1, 2 generators)
# W:=CoxeterGroup("E",6);   # (6 generators)
# W:=CoxeterGroup("E",7);   # (7 generators)
# W:=CoxeterGroup("E",8);   # (8 generators)
# W:=CoxeterGroup("F",4);   # (4 generators)
# W:=CoxeterGroup("H",3);   # (3 generators)
# W:=CoxeterGroup("H",4);   # (4 generators)

# b=B(<array>); / must be a positive element of "B", where "<array>" is an array (list) in the letters 1,2,...,k (k is number of generators).

#--------------------------------------------------------------------------------------------------------------------------------------------

LAlphaI:=function(b,I) # to compute the "<I>"-head of the element "b", where "I" is a subset of generators (array).
 	 local M,res,i,s;
	 M:=b.monoid; res:=M.Elt([]); i:=1;
	 while i<=Length(I) do
	       if M.IsLeftDescending(GarsideAlpha(b),I[i]) then s:=M.Elt([M.atoms[I[i]]]); res:=res*s; b:=s^-1*b; i:=1; else i:=i+1; fi;
  	 od;
  	 return [res,b];
	 end;

#--------------------------------------------------------------------------------------------------------------------------------------------

RAlphaI:=function(b,I) # to compute the "<I>"-tale of the element "b", where "I" is a subset of generators (array).
  	 return [ReversedWord(LAlphaI(ReversedWord(b),I)[2]),ReversedWord(LAlphaI(ReversedWord(b),I)[1])]; 
	 end;

#--------------------------------------------------------------------------------------------------------------------------------------------

AltDec:=function(b,J,I) # to compute the (J,I)-decomposition of "b", where "I,J" ae subsets of generators (arrays). 
	local h,M,i,dec;
	M:=b.monoid; h:=b; i:=1; dec:=[];
	while AsWord(h)<>[] do
	      if i=1 then Add(dec,RAlphaI(h,I)[2]); h:=h*RAlphaI(h,I)[2]^-1; i:=2;
	      else Add(dec,RAlphaI(h,J)[2]); h:=h*RAlphaI(h,J)[2]^-1; i:=1; fi;
	od;
	return Reversed(dec);
	end;

#--------------------------------------------------------------------------------------------------------------------------------------------

