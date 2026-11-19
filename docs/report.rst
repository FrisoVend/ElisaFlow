Rapport
=================

Bij formuleberekeningen staat onder **Assaygegevens** het resultaatlabel met
de berekening tussen haakjes, bijvoorbeeld **Berekening: Ratio (Monster / Cutoff)**.
**Monster** staat voor de extinctie van het betreffende monster.

.. _rapport-controleren:

Teststatus en waarschuwingen
----------------------------

**Geldig**, **Ongeldig** en **Niet bepaald** volgen uit de
:doc:`QC-beoordeling <validation>`.

Onder **Kwaliteitscontrole** staan per criterium de grens, beoordeelde waarde,
eenheid, CV (indien beschikbaar) en status. **Ext** betekent vóór blancoaftrek;
**ΔExt** betekent dat de blanco is afgetrokken. De waarde en CV volgen de
extinctieschaal van het betreffende criterium. Beide schalen omvatten eventuele
referentiecorrectie en verzadigingsomrekening. Zonder blanco staat er **Ext**.

**Rood** betekent dat een ingestelde PDF-waarschuwingsgrens is overschreden.
Bij **CV waarschuwing (%) = 20** worden CV's boven 20% rood; bij
**Terugvinding afwijking (%) = 20** geldt dit voor terugvindingen onder 80%
of boven 120%. Instellen kan onder :ref:`Rapport <instellingen-rapport>`.

Extinctiedata en markeringen
-----------------------------
Onder de metadata en QC-tabel staan de rauwe extincties, bijvoorbeeld **405 - 620 nm**.
Er zijn een aantal gevallen waar andere informatie of tekens in de tabel staan, deze zijn:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Markering
     - Betekenis
   * - **SAT**
     - Verzadigde meting zonder bruikbare omrekening.
   * - **UIT**
     - Well uitgeschakeld in SkanIt.
   * - **REF**
     - Ontbrekende of onbruikbare referentiemeting.
   * - **\***
     - Extinctie berekend uit een andere golflengte wegens verzadiging.
   * - **-**
     - Geen bepaalbare waarde of statistiek.

De tabel **Uitgesloten replicaten** geeft aan welke wells niet meetellen en waarom.
Bij uitsluiting door een replicaatfilter vermeldt de reden ook de gekozen
extinctieschaal en het oorspronkelijke replicaatgemiddelde.

Kalibratie en meettijdstip
--------------------------

Een gefitte curve toont de grafiek en de terugvinding van standaarden.
Ook vaste curven krijgen een grafiek, zonder fitstatistiek of terugvindingstabel.
Standaarden zijn daarvoor niet vereist. Het model en de coëfficiënten staan
in het legendavak van de grafiek. Een ingestelde **Voorbewerking** staat in
hetzelfde vak, in dezelfde tekstnotatie als **Berekening** bij de assaygegevens,
met de namen van **Monster** en de gebruikte controles. De y-as geeft de
extinctieschaal aan.
Mislukte passingen blijven vermeld, ook als een andere golflengte resultaten levert.

Patiëntresultaten
-----------------

De extinctie en bijbehorende CV volgen de gebruikte
berekening: **Ext** voor ruwe extinctie, **ΔExt** na blancoaftrek.
Bij kalibratievoorbewerking blijven deze korte labels staan; de vergelijking
van de voorbewerking staat in het legendavak van de bijbehorende kalibratiegrafiek.

De gerapporteerde concentratie en meetgrens zijn al gecorrigeerd voor 
verdunningen.

Een verzadigde meting kan een :ref:`concentratiegrens <rapportagebereik>`
opleveren, bijvoorbeeld **>300.0**, terwijl de plaat **SAT** blijft tonen.
