Program konsolowy do symulacji nieustalonej wymiany ciepła metodą elementów skończonych.\
Sposób użycia: 

`messim [-ip 2|3|4] [-o <outFileName> [-p]] <gridFile>`

Opcja `-ip` (lub `--integration-points`) decyduje o ilości punktów całkowania (2, 3 lub 4). \
Opcja `-o` (lub `--output`) decyduje o zapisie do pliku - domyślnie wynik symulacji wyświetlany jest w
terminalu, w postaci najmniejszej i największej temperatury ze wszystkich wierzchołków w danym kroku czasowym. \ 
Opcja `-p` (lub `--paraview-output`) decyduje o zapisie w formacie kompatybilnym z oprogramowaniem
ParaView, co daje możliwość głębszej wizualizacji procesu.
