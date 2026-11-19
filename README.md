# ElisaFlow

Windows-desktopapplicatie voor SkanIt ELISA-metingen, assayberekeningen,
kwaliteitscontrole en PDF-rapporten, vereist SkanIt 8.0 sessiebestanden.

**Disclaimer:** Deze software is bedoeld voor onderzoeksgebruik. Resultaten dienen
onafhankelijk te worden gevalideerd voordat ze worden gebruikt voor klinische,
diagnostische of andere beslissingen.

Open **ElisaFlow.exe** vanuit de volledige applicatiemap. Gebruik **Start meting**
voor metingen en **Assay beheer** voor instellingen.
Met **Heranalyse** kies je een opgeslagen experimentmap om opnieuw te analyseren.

## Documentatie

De Nederlandstalige Sphinx-documentatie is de centrale handleiding:

- [Gebruikershandleiding](docs/usage.rst): metingen en heranalyse.
- [Assay beheer](docs/assays.rst): protocollen en instellingen per scherm en tabblad.
- [Berekeningen](docs/calculation.rst): correcties voor rauwe data en formules.
- [Standaardcurven](docs/calibration.rst): Modellen en opties voor curve-fitten.
- [Kwaliteitscontrole en replicaatfilters](docs/validation.rst): QC-criteria en validatie.
- [Rapport](docs/report.rst): uitslagen, markeringen en teststatus.
- [Verwerking](docs/architecture.rst): volgorde en opslag.
- [API-referentie](docs/api.rst): Python, configuratie, command line en build-instructies.

Na het bouwen staat de handleiding in [build/docs/index.html](build/docs/index.html).

## Licentie

Deze software is beschikbaar onder de [MIT-licentie](LICENSE).
Copyright (c) 2026 Reinier Diagnostiek.
