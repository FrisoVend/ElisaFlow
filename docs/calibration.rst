Standaardcurven
===============

**Instellen:** :ref:`Berekening <instellingen-berekening>`.
Gebruik altijd het model, de weging en het bereik die voor de assay zijn vastgelegd.
Vaak geeft de fabrikant dit echter niet aan en alleen dat de fit goed moet zijn. Daarom hier
een korte uitleg van de belangrijkste factoren om deze keuze te maken.

Modellen
--------

Hieronder is :math:`x` de concentratie en :math:`y` de extinctie-invoer na eventuele
voorbewerking (bij vaste curves).

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Curvemodel
     - Parameters
     - Minimaal aantal verschillende standaardniveaus bij passing
   * - **4pl**
     - A, B, C, D
     - 4
   * - **5pl**
     - A, B, C, D, E
     - 5
   * - **Kwadratisch**
     - A, B, C
     - 3
   * - **Kubisch**
     - A, B, C, D
     - 4
   * - **Machtsfunctie met offset**
     - A, B, C
     - 3

.. math::

   y_{\mathrm{4PL}} = D + \frac{A-D}{1+(x/C)^B}
   \qquad
   y_{\mathrm{5PL}} = D + \frac{A-D}{(1+(x/C)^B)^E}

:math:`A` en :math:`D` zijn de lage en hoge asymptoot; :math:`B` bepaalt de richting en steilheid.
Bij :math:`E=1` wordt 5PL gelijk aan 4PL en ligt :math:`C` halverwege de respons.

.. caution::
   Er worden verschillende definities voor een 4PL en 5PL model gebruikt, wiskunding zijn ze allemaal identiek
   maar de waarde van de parameters is verschillend. De Q-koorst IgG assay gebruikt bijvoorbeeld de vorm
   :math:`y = A' + \frac{D' - A'}{1 + \exp{(B'(C' - \ln(x))})}`. 
   De waarden van de Q-koorts IgG kit moeten in dit geval als volgt omgerekend worden:

   * :math:`A = A'`
   * :math:`B = -B'`
   * :math:`C = e^{C'}`
   * :math:`D = D'`


.. math::

   y_{\mathrm{kwadratisch}} = A + Bx + Cx^2
   \qquad
   y_{\mathrm{kubisch}} = A + Bx + Cx^2 + Dx^3

Gefitte polynomiale curven zijn alleen stijgend of alleen dalend over het hele bereik (monotoon).

De **Machtsfunctie met offset** gebruikt:

.. math::

   y_{\mathrm{macht}} = C + A x^B

:math:`C` is de offset, :math:`A` de schaal en :math:`B` de exponent.
:math:`A` en :math:`B` mogen positief of negatief zijn, maar niet nul.
De curve is stijgend bij :math:`AB>0` en dalend bij :math:`AB<0`.
Bij een positieve exponent is de respons bij concentratie nul gelijk aan
:math:`C`, bij een negatieve exponent zijn alleen positieve concentraties
toegestaan. De machtsfunctie heeft geen verzadigingsplateau en kan daarom alleen
gebruikt worden voor standaarden rondom het stijgende en lineaire gedeelte.

Responsschaal en weging
-----------------------

* **Responsschaal:** kies **Vóór blancoaftrek** of **Na blancoaftrek** voor
  zowel curvepassing als vaste curven. Beide schalen zijn na eventuele
  referentiecorrectie en verzadigingsomrekening. Zonder blanco zijn ze gelijk.
* **Curve fitten uit standaarden:** standaard is **Vóór blancoaftrek**.
  De gekozen schaal geldt voor alle standaarden, controles en patiënten,
  ook bij het terugrekenen van concentraties.
* **Vaste parameters gebruiken:** kies de schaal waarop de parameters zijn
  bepaald. Wisselen tussen passen en vaste parameters behoudt de gekozen schaal;
  de parameters worden niet automatisch naar een andere schaal omgerekend.

De weging wordt bepaald door de **Wegingsexponent (1/y^p)** instelling, de fitter
minimaliseert de gewogen som van kwadratische afwijkingen:

.. math::

   w_i=y_i^{-p},
   \qquad
   \hat{\boldsymbol{\theta}} =
    \underset{\boldsymbol{\theta}}{\operatorname{arg\,min}}
    \sum_{i=1}^{n}w_i
    \left(y_i - f(x_i;\boldsymbol{\theta})\right)^2

:math:`p=0` is ongewogen, :math:`p=1` geeft :math:`1/y` en :math:`p=2` geeft :math:`1/y^2`.
Hier is :math:`y` de individuele standaardrespons op de gekozen extinctie schaal.
Een andere responsschaal verandert dus ook de gewichten.
Ook decimale en negatieve exponenten zijn mogelijk. Bij ELISA assays is :math:`1/y^2` het meest voorkomend, maar dit moet
per assay getest worden door de gemiddelde extinctie (Gem) tegen de standaard-deviatie (SD) te plotten voor alle standaarden
en te fitten met :math:`\log(\mathrm{SD}) = a + b \log(\mathrm{Gem})` (na referentiegolflengte aftrekking maar voor blancoaftrekking), 
waarna :math:`1/y^{2b}` het beste resultaat geeft.
Vaste parameters gebruiken geen weging of standaarden voor de berekening omdat de curve al bekend is.

Voorbewerking bij een vaste standaardcurve
------------------------------------------

**Instellen:** :ref:`Berekening → Voorbewerking <instellingen-berekening>`.

Bij vaste curves die aangeleverd worden door de kit-fabrikant is het mogelijk om een optionele formule
in te vullen, die extinctie metingen omrekent voordat de concentratie wordt afgelezen. Deze formule geldt voor
alle controles en patiënten op die golflengte. Als dit nodig is wordt het altijd aangegeven in de kit.

**Voorbewerking:** is alleen beschikbaar bij vaste parameters en wordt
toegepast op de gekozen schaal vóór concentratieberekening.

.. note::
   Een voorbeeld van voorbewerking bij een vaste curve is het Q-koorts IgG assay.
   Waar dag-tot-dag variaties worden gecorrigeerd met de F-factor,
   (de RV-waarde uit de kit gedeeld door de gemiddelde waarde van het Standaard monster).

Terugvinding
------------

Voor een standaardconcentratie geldt:

.. math::

   \mathrm{Terugvinding}(\%) =
   100\,\frac{\overline{c}_{\mathrm{terugberekend}}}{c_{\mathrm{nominaal}}}

.. _rapportagebereik:

Bereik en verdunning
--------------------

Een gefitte curve rapporteert binnen het interval van de behouden
standaarden, beperkt tot het ingestelde **Concentratiebereik**. Buiten dit interval wordt niet geëxtrapoleerd. Dit geldt bij één
en bij twee analysegolflengten. Een ingesteld bereik kan dit rapportagebereik
verkleinen, maar niet vergroten. 

Een verzadigde patiëntmeting zonder bruikbare omrekening kan eveneens een
grens opleveren. Verzadiging ligt boven de betrouwbare ruwe meetsignalen op
dezelfde golflengte. Het hoogste signaal van de behouden, onverzadigde
standaardmetingen dient daarom als ondergrens voor het ruwe patiëntsignaal.
Deze signaalgrens wordt met de eigen referentiemeting van de patiënt en de
responsschaal van de curve gecorrigeerd. Bij een stijgende curve volgt
**>bovengrens** als de gecorrigeerde signaalgrens ten minste de respons bij de
bovengrens bereikt; bij een dalende curve kan zo **<ondergrens** volgen.
De concentratiegrens wordt niet als numerieke uitslag of gemiddelde gebruikt.

Voorbewerking van de vorm :math:`aE+b` wordt meegenomen, inclusief een
negatieve helling en controleverwijzingen. Zonder betrouwbare
standaardmetingen, bij andere vormen van voorbewerking, bij een onbruikbare
referentie of als de signaalgrens het curve-eindpunt niet ondersteunt, blijft
de uitslag onbeschikbaar. Ook een vaste curve vereist deze standaardmetingen
en een ingesteld rapportagebereik om verzadiging als grens te rapporteren.
Een geslaagde omrekening of een numerieke uitslag op de andere golflengte
houdt voorrang. Bij strijdige replicaten, bijvoorbeeld een numerieke uitslag
en een verzadigde meting boven het bereik, wordt geen gemiddelde gerapporteerd.

Voor patiënten wordt na inversie en bereikcontrole met de SKAX-verdunningsfactor
:math:`d` vermenigvuldigd als

.. math::

   c_{\mathrm{rapport},i}=d\,c_i.

Bij een verdunning 1:10 is :math:`d=10`; zonder ingevulde factor geldt 1 en verschillende patiëntverdunningen blijven apart.

.. _golflengteselectie:

Twee analysegolflengten
-----------------------

**Instellen:** :ref:`Meting <instellingen-meting>` en het bijbehorende
golflengtetabblad onder **Berekening**.

Elke golflengte kan een eigen vaste of gefitte curve gebruiken.

Onder de **Omschakelconcentratie (vóór verdunningscorrectie)** heeft de
voorkeurscurve voorrang. Zodra één behouden replicaat die grens bereikt of
overschrijdt, gaat het volledige monster naar de andere curve. 

Een onbruikbare voorkeurscurve, verzadiging of een resultaat buiten
het bereik kan ook tot gebruikt van de andere curve leiden. Als overschrijding van de
omschakelconcentratie de andere curve vereist en die onbruikbaar is, blijft het
resultaat onbeschikbaar.
