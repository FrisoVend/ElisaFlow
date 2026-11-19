API-referentie
==============

Overzicht voor installatie, beheer en gebruik vanuit Python.
Voor de dagelijkse bediening zie :doc:`usage`; de volgorde van inlezen,
berekenen, kwaliteitscontrole en rapportage staat in :doc:`architecture`.

Installeren en starten
----------------------

Gebruik Python 3.15.2 of nieuwer. Voer vanuit de projectmap uit:

.. code-block:: powershell

   py -3.15 -m venv .venv
   .\.venv\Scripts\python.exe -m pip install -e ".[ui,dev,docs]"
   .\.venv\Scripts\python.exe -m elisaflow

.. _applicatiepaden:

Paden en instellingen
---------------------

``app.toml`` staat bij de executable, of bij broncodegebruik in de werkmap.
Relatieve paden worden vanaf dat bestand bepaald:

.. code-block:: toml

   [paths]
   protocols = "protocols"
   configs = "configs"
   results = "resultaten"

   [access_control]
   password_file = ".secrets/assay-editor.hash"

De protocolmap bevat de SkanIt-sjablonen; de configuratiemap bevat de
bijbehorende assay-instellingen. Wijzig deze via :doc:`assays`.

.. _bestandsbeheer:

Bestanden beheren
-----------------

Sluit ElisaFlow vóór handmatige bestandswijzigingen.

* **Assay archiveren:** verplaats de bijbehorende JSON uit ``configs/`` en het
  SKAX-protocol uit ``protocols/`` naar een archiefmap. Controleer welke
  protocollen de JSON bevat. Verwijder zo nodig de Windows-eigenschap
  **Alleen-lezen**.
* **Wachtwoord herstellen:** verwijder het ingestelde wachtwoordbestand en stel
  na herstart een nieuw wachtwoord in via **Wachtwoord beheren**.
* **Storingen onderzoeken:** het logbestand staat standaard in
  ``%LOCALAPPDATA%\ElisaFlow\ElisaFlow.log``.

.. _command_line:

Command line
-------------

Vanuit de gebouwde applicatiemap kun je een meting of heranalyse starten:

.. code-block:: powershell

   .\ElisaFlow.exe run --protocol ".\protocols\Hib_IgG.skax"
   .\ElisaFlow.exe run --experiment "C:\ElisaFlow\resultaten\vorige-meting"

Zonder opdracht opent de desktop; ``settings`` opent **Assay beheer**.
Iedere meting of heranalyse vraagt de analistnaam. Alle opties staan bij
``.\ElisaFlow.exe run --help``.

Voor heranalyse kies je een map met één voltooide SKAX-sessie en één
``*_settings.json`` met de gebruikte assay-instellingen.

Vanuit de broncode gebruik je dezelfde opdrachten en argumenten na
``.\.venv\Scripts\python.exe -m elisaflow``.

Gebruik vanuit Python
---------------------

``run_workflow`` voert het hele proces uit: invoer bewaren, meetgegevens
inlezen, berekenen, kwaliteitscontrole uitvoeren en een PDF maken.
Bij een live meting opent het ook SkanIt en wacht het op het uitvoeren van de meting.
Dit voorbeeld heranalyseert een bestaande meting:

.. code-block:: python

   from pathlib import Path

   from elisaflow.configuration.application import ApplicationSettings
   from elisaflow.runs.workflow import run_workflow

   config = ApplicationSettings.from_file("app.toml")
   report = run_workflow(
       config,
       analyst="Naam analist",
       experiment_dir=Path("resultaten/vorige-meting"),
   )

``report`` is het pad naar de nieuwe PDF. Gebruik voor een live meting
``protocol_template=Path("protocols/Hib_IgG.skax")`` in plaats van
``experiment_dir``.

Voor alleen de analyse gebruik je ``AssayAnalyzer`` uit
``elisaflow.analysis.assay``: ``result = AssayAnalyzer(data, settings).analyze()``.
Hierbij zijn ``data`` en ``settings`` al ingelezen meetgegevens en
assay-instellingen. ``result.samples`` bevat de geselecteerde monsters,
``result.measurement_results`` de uitkomsten per golflengte en
``result.validation_results`` de QC-beoordelingen van de geselecteerde resultaten.

Bouwen en controleren
---------------------

Maak de Windows-applicatie:

.. code-block:: powershell

   powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\build_windows.ps1

Kopieer de volledige map ``dist/build_<timestamp>/ElisaFlow/`` naar een
schrijfbare locatie op de doelcomputer. Python en de offline handleiding
worden meegeleverd. Voor live metingen is SkanIt nodig. Controleer daar een volledige meting en heranalyse met bekende resultaten.

Testen en documentatie:

.. code-block:: powershell

   $env:QT_QPA_PLATFORM = "offscreen"
   $env:MPLBACKEND = "Agg"
   .\.venv\Scripts\python.exe -m ruff check src tests scripts
   .\.venv\Scripts\python.exe -B -m unittest discover -s tests -p "test_*.py" -q
   .\.venv\Scripts\python.exe -m sphinx -a -W -b html docs build/docs
