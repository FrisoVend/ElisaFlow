Berekeningen
============

De instelling **Resultaatmethode** bepaalt of een formule of een standaardcurve
het resultaat berekent. Iedere behouden well wordt afzonderlijk berekend,
daarna wordt het gemiddelde per sample berekend. De
:doc:`verwerkingsvolgorde <architecture>` ligt vast.

.. _correcties:

Referentie, verzadiging en blanco
---------------------------------

**Instellen:** :ref:`Meting <instellingen-meting>`.

Een referentiegolflengte (bijna altijd 620nm) wordt voor iedere well eerst afgetrokken:

.. math::

   \mathrm{Ext}_{\lambda,\mathrm{ref}} = \mathrm{Ext}_\lambda - \mathrm{Ext}_{\mathrm{referentie}}

Bij **Berekenen uit een andere golflengte** geldt voor een verzadigde meting:

.. math::

   \mathrm{Ext}_{\mathrm{doel}} = b + a \mathrm{Ext}_{\mathrm{bron}}

Hier zijn :math:`a` de **Helling** en :math:`b` de **Intercept**.
De omrekening gebruikt een onverzadigde bronmeting van een andere golflengte. 
Na aftrekking van de referentiegolflengte maar vóór blanco-aftrekking.

Zonder bruikbare omrekening blijft een verzadigde meting uitgesloten als controle of standaard. 
Bij patiënten kan een verzadigde meting wel een :ref:`concentratiegrens <rapportagebereik>` opleveren, 
omdat een verzadigde extinctie altijd hoger is dan de hoogste onverzadigde standaard op die golflengte.
De onbetrouwbare extinctiewaarde zelf wordt daarbij niet gebruikt, en de plaatweergave houdt
de markering **SAT**.

Als blanco's zijn opgenomen, wordt per analysegolflengte het
gemiddelde afgetrokken:

.. math::

   \Delta \mathrm{Ext}_i = \mathrm{Ext}_i - \overline{\mathrm{Ext}}_{\mathrm{blanco}}

De waarden vóór en ná blancoaftrek blijven beschikbaar. Voor extinctie-QC,
curve-invoer en formule-invoer wordt de afzonderlijk gekozen of de extincite vóór of ná 
blanco-aftrekking gebruikt worden.

.. _resultaatformules:

Resultaatformules
-----------------

**Instellen:** :ref:`Berekening → Patiënten → Formule <instellingen-berekening>`.

Kies bij iedere formule een **Responsschaal**: **Vóór blancoaftrek** of
**Na blancoaftrek**. Patiënten en iedere afwijkende controleformule hebben een
eigen keuze. De eigen meetwaarde en alle controles in een formule
gebruiken de schaal van die formule. Met een blanco in het protocol is **Na blancoaftrek** de standaard,
zonder blanco is dit **Vóór blancoaftrek** en zijn beide schalen gelijk. Iedere keuze wordt
expliciet opgeslagen.

Met **Herstellen** bij de controleformule gebruik je opnieuw de patiëntformule,
terwijl de eigen invoerschaal behouden blijft.

Een voorbeeld van een berekening is de AP50-formule, die de extinctie van de negatieve controle op 0% zet en de positive controle op 100%:
In het formulevenster is dit:

.. code-block:: text

   100 * (Patient - Negatief) / (Positief - Negatief)

Formule-assays vereisen een patiëntverdunningsfactor van 1, omdat er geen linearie relatie 
is tussen de extincite en de uitslag. 



