# Thema 6: Anwendungsbeispiele

... wird noch vervollständigt ...

[TOC]


<!------------------------------------------------------------------------------
Bewährte Programmierpraktiken
------------------------------------------------------------------------------->
## Bewährte Programmierpraktiken

In dem [Begleitmaterial](Begleitmaterial) befinden sich bereits Beispiele für bewährte Programmierpraktiken. Denn das Einhalten von Konvention macht die Entwicklung von Programmcode, für alle Beteiligten, sehr viel leichter. Genauso wie das wissenschaftliche Arbeiten zum Verfassen von fundierten Berichten gehört, so ist auch der Programmierstil für die Entwicklung von Programmcode maßgeblich. Dazu folgende Anmerkungen.

1. **Dokumentation:** Der Programmcode lässt sich durch Kommentare in Sinnabschnitte unterteilen. Außerdem bietet es sich an, mit einem sog. Docstring, kurze Funktionsbeschreibungen vorzunehmen.

2. **Vermeidung von Redundanz:** Werte die mehrmals vorkommen, sollten als Variablen verwendet werden. Routinen welche mehrmals vorkommen, sollten als Funktionen ausgelagert werden.

3. **Verwendung von Parametern:** Bei einer Simulation gibt es immer Parameter. Es bietet sich dafür an, die Hauptfunktion mit sog. Keyword-Arguments zu versehen, um dadurch Voreinstellungen für die möglichen Parameter zu treffen.

4. **Typisierung:** Durch das sog. Type-Hinting kann bei Programmiersprachen mit impliziter Typisierung explizit auf den Datentyp hingewiesen werden. Dies ist gerade bei Funktionsargumenten und Rückgabewerten sehr nützlich, da so Missverständnisse vermieden werden.

5. **Handhabung von großen Datenmengen:** Bei Skriptsprachen wie Matlab oder Python erfolgt der Funktionsaufruf meistens durch sog. Wertparameter (call by value) und nicht durch sog. Referenzparameter (call by reference), wodurch das gesamte Funktionsargument bei dem Funktionsaufruf kopiert wird. Da Speicheroperationen aber in der Regel sehr langsam sind, ist es in diesem Fall günstiger, große Objekte – wie z. B. die Differenzenmatrizen – global zu definieren, um sie nicht ständig neu zu initialisieren.


<!------------------------------------------------------------------------------
Taylor-Green-Wirbel
------------------------------------------------------------------------------->
## Taylor-Green-Wirbel

...

### Analytische Lösung

### Visualisierungsbeispiel

### Untersuchung
