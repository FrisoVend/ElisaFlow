Verwerking
=========================

Volgorde
--------

1. **Invoer controleren.** De plaatindeling en metingen worden uit de
   opgeslagen SKAX gelezen. Blanco-, controle- en standaardwells moeten bij
   elke protocolgolflengte zijn opgenomen. Elke meting heeft een fysieke
   wellcoördinaat en een golflengte.
2. **Kanalen voorbereiden.** Referentiemetingen worden per well afgetrokken.
   Eventuele verzadigingsomrekening gebruikt de rauwe extinctiemeting,
   na referentieaftrek. Lege patiëntposities worden bepaald vóór blancoaftrek.
3. **Per analysegolflengte rekenen.** Onbruikbare wells worden uitgesloten;
   daarna volgen blancoverwerking, replicaatfiltering en de formule of
   kalibratie. Optionele kalibratievoorbewerking wordt per curve gedaan.
   Berekende waarden bewaren de gebruikte extinctieschaal. Gefitte
   curven rapporteren alleen binnen het bereik van de behouden standaarden,
   eventueel verder beperkt door het ingestelde concentratiebereik.
4. **Monsterresultaten kiezen.** Bij twee analysegolflengten wordt per monster
   één golflengte gekozen. Verdunningscorrectie hoort bij het berekende
   patiëntresultaat, de omschakelgrens wordt vóór die correctie toegepast.
5. **QC en rapport maken.** Criteria per golflengte en voor geselecteerde
   resultaten bepalen samen de teststatus. De PDF toont resultaten,
   waarschuwingen en uitsluitingen.

De rekenregels staan in :doc:`calculation`, en uitleg over het rapport in
:doc:`report`.

Uitslag
---------------------

Een fout die de berekening onmogelijk maakt, stopt de verwerking.
Een QC-afkeuring kan wel een rapport met berekende waarden opleveren, wel wordt
aangegeven waar de fout zit. Controleer daarom het rapport handmatig voordat je resultaten overneemt.
