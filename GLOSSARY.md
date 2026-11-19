# ElisaFlow

ElisaFlow verwerkt ELISA-metingen uit SkanIt tot berekende uitslagen,
kwaliteitsbeoordelingen en analyserapporten.

## Taal

**Patiëntpagina**:
Het deel van het analyserapport met patiëntuitslagen, bijbehorende meetwaarden
en variatiecoëfficiënten.

**Patiëntresultaat**:
De uitslag voor een patiëntmonster op basis van de behouden replicaten en de
uiteindelijk gekozen analysegolflengte, met de toepasselijke verdunningscorrectie.

**Resultaatkwalificatie**:
De betekenis van een uitslag als numeriek resultaat binnen het rapportagebereik,
als onder- of bovengrens, of als niet beschikbare uitslag. Een grens is geen
numerieke concentratie en draagt niet bij aan een resultaatgemiddelde of -CV.

**Extinctieschaal**:
De betekenis van extinctiewaarden volgens de toegepaste referentieaftrek,
verzadigingsomrekening, blancoverwerking en eventuele kalibratievoorbewerking.

**Ruwe extinctie (Ext)**:
De extinctie vóór blancoaftrek, na eventuele referentieaftrek en
verzadigingsomrekening.

**Blanco-gecorrigeerde extinctie (ΔExt)**:
De extinctie na aftrek van het bruikbare blancogemiddelde, ook wanneer dat
gemiddelde nul is.

**Curve-invoer**:
De extinctiewaarde per well die de kalibratie gebruikt, op de ingestelde
responsschaal en na eventuele voorbewerking.

**Rapportagebereik**:
Het concentratie-interval waarbinnen een kalibratie een numerieke uitslag
rapporteert. Bij een gefitte curve is dit het bereik van de behouden positieve
standaarden, eventueel verder beperkt door het ingestelde concentratiebereik;
bij vaste parameters geldt het ingestelde bereik zonder vereiste standaarden.

**Fitdomein**:
Het concentratie-interval van alle behouden standaarden, inclusief eventuele
nulstandaarden, waarover een polynomiale curvepassing monotoon moet zijn.

**Monotone inversietak**:
Het monotone deel van een polynomiale curve waarop een extinctiewaarde
eenduidig naar een concentratie wordt omgerekend.

**Extinctie-QC**:
De kwaliteitsbeoordeling van extinctiewaarden op de daarvoor geldende schaal;
die kan verschillen van de curve-invoer van het patiëntresultaat.

**Validatiecriterium**:
Een onder- en/of bovengrens voor een bepaalde grootheid van een blanco,
benoemde controle of standaardconcentratie. Bij extinctie hoort daarbij een
expliciete keuze vóór of na blancoaftrek.

**Experimentmap**:
Een map met de bewaarde SkanIt-sessie en instellingensnapshot van een meting
of heranalyse, eventueel aangevuld met een analyserapport.

**Instellingensnapshot**:
De vastgelegde assay-instellingen die voor een meting of heranalyse zijn
toegepast, onafhankelijk van latere wijzigingen in de assaybibliotheek.
