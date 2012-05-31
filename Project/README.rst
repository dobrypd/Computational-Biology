Main Project
----
----

Requirements
----
Blast+
can be downloaded from here:
http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Sequences:
new organism:
http://tofesi.mimuw.edu.pl/~bartek/WBO/yeast_relative.fa
S. cerevisiae sequences:
http://tofesi.mimuw.edu.pl/~bartek/WBO/yeast.fa
positions:
http://tofesi.mimuw.edu.pl/~bartek/WBO/yeast.bed
accotations:
http://tofesi.mimuw.edu.pl/~bartek/WBO/yeast_annotations.go

Task
----
link: http://tofesi.mimuw.edu.pl/wbo/?p=196

To zadanie nieco różni się od poprzednich. O ile do tej pory starałem się dość dokładnie specyfikować jakie ma być wejście, wyjście i parametry funkcji, które Państwo pisali, o tyle projekt zaliczeniowy wymaga od Państwa podejmowania w takich sprawach decyzji opartych o to czego Państwo się do tej pory nauczyli i ew. własne doświadczenia w trakcie rozwiązywania zadania. Punktowanie tego projektu też będzie wyglądało nieco inaczej: po wysłaniu gotowego rozwiązania będą Państwo spotykać się ze mną (na zajęciach albo konsultacjach) i opowiadać co udało się zrobić, jak zostało to zrobione i dlaczego właśnie tak a nie inaczej. W razie gdyby pojawiły się wątpliwości, oczywiście jestem do Państwa dyspozycji i chętnie będę udzielał wszelkich informacji pomocnych w zrozumieniu zadania.
Opis funkcjonalny nowej sekwencji genomowej
Naszym zadaniem jest opisanie pewnych aspektów funkcjonalnych nowerj sekwencji genomu drożdża, odległego krewnego drożdży piekarskich (S. cerevisiae) znanego organizmu modelowego. Mamy do dyspozycji samą sekwencję DNA nowego organizmu a także zestaw plików opisujących znany już genome S. cerevisiae :
* sekwencja genomowa,
* opis położenia genów,
* przypisanie genów do funkcji wg. Gene Ontology.
Nasze zadanie podzielimy na 5 elementów, każdy z nich jest oceniony niezależnie, ale jak łątwo zauważyć ich sekwencja jest nieprzypadkowa. Nie jest obowiązkowe wykonanie wszystkich zadań, ani też nie jest konieczne ścisłe trzymanie się kolejności, ale niektóre zadania będą wymagały wyników zadań wcześniejszych.
Wyszukanie prawdopodobnych genów ortologicznych (30 pkt) 
    Na podstawie genomu  i pozycji genów S. cerevisiae, wyszukaj prawdopodobne lokalizacje (początek, koniec i nić) genów w nowym genomie. Wynik powinien zawierać przypisania homologów i jakąś statystyczną miarę istotności każdej predykcji
Określenie prawdopodobnych funkcji genów na podstawie PFAM (20 pkt) 
    Każdy z przewidzianych genów może zostać przeanalizowany na obecność domen białkowych wg. bazy PFAM. Chcemy wiedzieć które geny mają przewidziane domeny PFAM i czy te predykcje są zgodne z predykcjami domen w ich ortologach w S. cerevisiae.
Wyszukanie potencjalnych obszarów regulatorowych w nowym genomie (10 pkt).
    Zakładamy, że każdy gen jest regulowany transkrypcyjnie przez czynniki transkrypcyjne wiążące się do obszaru DNA rozciągającego się od początku transkrypcji w kierunku 5′ aż do następnego genu (obszar międzygenowy “upstream”). Znajdź pozycje takich obszarów regulatorowych w nowym genomie i określ pary obszarów przypisane do wszystkich par “ortologów”.
Wyszukiwanie zachowanych ewolucyjnie miejsc wiązania czynników transkrypcyjnych (30 pkt)
Korzystając z bazy dancych JASPAR core fungi określ występujące w obszarach regulatorowych wystąpienia miejsc wiązania czynników transkrypcyjnych, które są zachowane w ewolucji pomiędzy tymi dwoma gatunkami drożdży.
Statystycznie nadreprezentowane pary motyw-funkcja (20 pkt).
    Możemy zastanowić się, które motywy występują w obszarach regulatorowych genów o szczególnych funkcjach (wg. PFAM i Gene Ontology). Jeśli jakiś motyw występuje w promotorach genów o określonej funkcji częściej niż spodziewalibyśmy się tego w losowym przypadku, możemy podejrzewać jego rolę w regulacji tejże funkcji. Zaproponuj statystykę mierzącą “nadreprezentację” i podaj nadreprezentowane pary motyw-funkcja.
