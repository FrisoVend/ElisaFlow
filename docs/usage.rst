Gebruikershandleiding
=====================

Starten
-------

Open **ElisaFlow.exe** vanuit een folder met schrijfrechten. Voor metingen zijn SkanIt 8.0 en een PDF-lezer nodig.

Met **Open documentatie** rechtsonder, onder de tabbladen, open je deze
handleiding in de browser.

Stel een assay eerst in via :doc:`assays`, als er nog geen enkele geconfigureerd is.

.. _live-meting:

Een meting uitvoeren
--------------------

1. Klik op **Start meting**, selecteer de assay en vul je naam in.
2. Voer de meting uit in de SkanIt-sessie die ElisaFlow opent en wacht tot de PDF opent. Alle bestanden worden automatisch opgeslagen.

**Verwerking annuleren** stopt het wachten of de verdere verwerking.
Een lopende apparaatmeting moet je zelf in SkanIt stoppen.

.. _opgeslagen-meting:

Opslag en heranalyse
--------------------
Elke meting krijgt een eigen map onder ``resultaten/``, ingedeeld op jaar en week.
De map bevat de PDF, de voltooide SKAX-sessie en de gebruikte instellingen.

Bewaar de volledige map. Een heranalyse gebruikt de opgeslagen invoer en
instellingen met de analysecode van de huidige ElisaFlow-versie. De uitvoer daarvan
komt in een nieuwe map en de oorspronkelijke meting blijft behouden.

1. Klik onder **Start meting** in de sectie **Opgeslagen meting** op **Heranalyse**.
2. Kies de experimentmap met de voltooide SKAX-sessie en de opgeslagen instellingen.
3. Vul je naam in en controleer het nieuwe PDF-rapport.

Heranalyse is ook beschikbaar via de :ref:`command line <command_line>`.
