Kwaliteitscontrole en replicaatfilters
======================================

**Instellen:** :ref:`QC-criteria <instellingen-qc>`.
Een rode rapportwaarde is een waarschuwing; de ingestelde QC-grenzen bepalen
de acceptatie met onafgeronde waarden.

Welke waarde wordt beoordeeld?
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Criterium
     - Beoordeelde waarde
   * - **Blanco → Extinctie**
     - Het bruikbare blancogemiddelde vóór blancoaftrek, na eventuele referentiecorrectie.
   * - **Controle of standaard → Extinctie**
     - Extinctie op de bij het criterium gekozen schaal: vóór of na
       blancoaftrek, steeds vóór kalibratievoorbewerking.
   * - **Controle of standaard → Resultaat**
     - Het gemiddelde van de berekende replicaatresultaten.

Controles worden herkend aan hun exacte naam en standaarden aan hun nominale
concentratie. Per controle of standaard, en bij twee analysegolflengten per
golflengte, is de extinctieschaal afzonderlijk instelbaar. Dezelfde keuze geldt
voor minimum, maximum en de extinctie-CV in het rapport.

Bij nieuwe criteria is **Na blancoaftrek** de standaard als het protocol een
blanco bevat; anders is dit **Vóór blancoaftrek**. De gekozen schaal wordt
expliciet opgeslagen.

Bij twee analysegolflengten dragen alle ingestelde **Meting**-criteria bij aan
de totale beoordeling, ook als die golflengte geen patiëntresultaten levert.
**Resultaat** beoordeelt de geselecteerde uitslagen.

De variatiecoëfficiënt wordt berekend als :math:`100s/|\overline{x}|`, met
:math:`s` de steekproefstandaardafwijking. De CV is dus ook bij een negatief
gemiddelde niet-negatief. Bij een gemiddelde van nul is de CV niet bepaald.

Statussen en resultaten buiten bereik
-------------------------------------

* **OK:** waarde valt binnen alle vastgestelde grenzen.
* **Afgekeurd:** een grens wordt overschreden.
* **Niet bepaald:** de beoordeling is niet mogelijk, bijvoorbeeld door ontbrekende metingen

De totale teststatus is **Ongeldig** zodra één criterium faalt.
De status is **Niet bepaald** als een beoordeling ontbreekt
of geen QC-criteria zijn ingesteld. **Geldig** vereist dat alle ingestelde
criteria voldoen.

Lege wells
-------------------------

Patiëntwells met een bruikbare extinctie onder 0,009, vóór blancoaftrek, gelden
als lege stripposities. Bij meerdere golflengten wordt een monstergroep alleen als
leeg weggelaten als alle metingen op beide kanalen onder die grens liggen.
Dit leegcriterium geldt niet voor blanco's, controles of standaarden.

.. _replicaatfilter:

Replicaatfilter
---------------

**Instellen:** :ref:`Monsterinstellingen <instellingen-monsters>`.

Kies bij ieder filter een **Responsschaal**: **Vóór blancoaftrek** of
**Na blancoaftrek**. Patiënten, iedere controle en iedere standaardconcentratie
hebben een eigen keuze, onafhankelijk van berekening en QC. Met een blanco in
het protocol is **Na blancoaftrek** de standaard; anders **Vóór blancoaftrek**.
De gekozen schaal wordt altijd opgeslagen.

Het filter werkt per golflengte op de gekozen schaal, vóór formules en
kalibratievoorbewerking. Per well geldt:

.. math::

   D_i(\%)=100\,\frac{|E_i-\overline{E}|}{|\overline{E}|}

Hier is :math:`E_i` de extinctie op de gekozen schaal. Beide schalen zijn na
eventuele referentiecorrectie en verzadigingsomrekening. Zonder blanco zijn ze gelijk.

Wells boven de ingestelde afwijkingsgrens van het gemiddelde worden 
uitgesloten. Standaard is maximaal één filteruitsluiting toegestaan.
Uitgesloten wells blijven zichtbaar in het rapport met hun reden.
