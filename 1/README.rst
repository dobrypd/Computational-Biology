First homework:
---------------

Zadanie domowe na dziś:
Napisać funkcję, która dostaje na wejściu plik z listą sekwencji nukleotydowych w formacie FASTA. Na wyjściu powinniśmy otrzymać listę sekwencji, w postaci obiektów SeqRecord, opisanych tymi samymi identyfikatorami co na wejściu. Każda sekwencja na wyjściu powinna hybrydyzować z odpowiednią sekwencją z wejścia i nie hybrydyzować z żadną inną sekwencją wejściową (Przypominam o komplementarności i odwracaniu). Algorytm może być zupełnie naiwny (brutalne wyliczenie, np. przy pomocy zbiorów), ale powinien rozpoznawać sytuację gdy nie istnieje rozwiązanie (i np. zgłaszać wyjątek).
Rozwiązanie pełne powinno rozważać 3 możliwości:
W przypadku podania parametru k, wszystkie sondy powinny mieć długość k.
k jest parametrem opcjonalnym. W przypadku braku podania k, funkcja zwraca rozwiązanie dla minimalnego k, dla którego rozwiązanie istnieje.
Zamiast k można podać też parametr f opisujący więzy na zwracane sondy. Ten parametr powinien być funkcją w pythonie, która dla dowolnej sekwencji zwraca, jako wartość logiczną, czy może ona wystąpić w rozwiązaniu, czy nie. Przykładowa funkcja może np. ograniczać dozwolone temperatury topnienia (def) dla sekwencji w rozwiązaniu. Można też rozważyć opcję, w której długość rozwiązań nie musi być identyczna.
Np. jeśli funkcja f jest zdefiniowana następująco: lambda(s):len(s)==5, to wynik powinien być identyczny do podania parametru k=5 (choć czas działania może być dłuższy…).
