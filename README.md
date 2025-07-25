# Simulation et Optimisation d'un Procédé de Chromatographie Préparative en Bioraffinerie

## 📜 Contexte du Projet

Ce projet s'inscrit dans le cadre de la valorisation de la biomasse lignocellulosique, une ressource clé pour la production de biocarburants et de molécules biosourcées. Après une étape d'hydrolyse acide pour extraire la cellulose, on obtient un hydrolysat contenant des molécules d'intérêt (sucres comme le glucose) mais aussi des sous-produits (acide acétique, sels).

L'objectif de ce projet est de modéliser et d'optimiser un procédé de **chromatographie préparative** pour séparer efficacement ces composés, afin de maximiser la valorisation de chaque fraction. L'étude se concentre sur l'utilisation d'une résine anionique forte, couramment employée dans l'industrie pour ce type de séparation.

***

## 🎯 Objectifs

Le projet se décompose en trois objectifs principaux :
1.  **Étude Expérimentale :** Réaliser des élutions en colonne pour caractériser le comportement de chaque composé (NaCl, glucose, acide acétique) à différents débits.
2.  **Modélisation & Estimation de Paramètres :** Utiliser un simulateur MATLAB pour modéliser les profils d'élution expérimentaux et estimer les paramètres clés du modèle (porosité du lit, nombre d'étages théoriques - NET).
3.  **Optimisation de Procédé :** Simuler un procédé industriel multi-colonnes **ISMB** (Improved Simulated Moving Bed) à l'aide d'un second programme MATLAB, et optimiser ses conditions opératoires pour une séparation binaire donnée, afin de maximiser la productivité tout en minimisant la consommation d'éluant.

***

## 🧑‍💻 Ma Contribution

Au sein du **Groupe N°1**, ma mission était centrée sur l'étude expérimentale du **NaCl** et l'optimisation de la séparation du mélange binaire **NaCl/Glucose**. Mon travail, détaillé dans le rapport `rapport_rattrapage_EI_Chloé_Lucas_Aymen_Awainia.pdf`, a couvert l'ensemble du processus :

### 1. Phase Expérimentale
* **Réalisation des essais :** J'ai mené les expériences d'élution simple pour le NaCl sur une colonne de 150 mL. J'ai effectué les manipulations pour quatre débits distincts : 4, 3, 2 et 1 BV/h (Bed Volume/heure).
* **Analyse des échantillons :** J'ai réalisé le suivi des fractions collectées en sortie de colonne par conductimétrie, après avoir établi une courbe d'étalonnage pour convertir les mesures de conductivité en concentration massique (g/L).
* **Traitement des données :** J'ai traité les données brutes, en effectuant les corrections nécessaires (notamment pour les perturbations observées à 1 BV/h) et en appliquant un bilan matière pour garantir la cohérence des profils d'élution.

### 2. Modélisation et Simulation (MATLAB)
* **Estimation des paramètres :** En utilisant le script `Simu_Chromato_Batch_EI.m`, j'ai ajusté les paramètres du modèle (porosité et NET) pour que les courbes simulées correspondent aux profils expérimentaux obtenus pour le NaCl à chaque débit. Les valeurs identifiées sont :
    * **Porosité ($ε$) :** 0.55 (constante pour tous les débits).
    * **NET :** 100 (à 4 BV/h), 125 (à 2 BV/h), 150 (à 3 BV/h).
* **Analyse de la rétention :** J'ai utilisé l'ensemble des paramètres (les miens pour le NaCl et ceux des autres groupes pour le glucose et l'acide acétique) pour simuler la séparation d'un mélange ternaire et comparer la rétention des différents composés, confirmant que le sel est non retenu, suivi du glucose, puis de l'acide acétique.

### 3. Optimisation du Procédé ISMB
* **Simulation ISMB :** À l'aide du script `Simu_ISMB_2F_4Z_EI.m`, j'ai mené l'optimisation de la séparation du mélange binaire **NaCl/Glucose** pour les quatre débits étudiés.
* **Détermination des conditions optimales :** Par une approche itérative, j'ai ajusté les volumes de la séquence d'injection (BV1, BV2, BV3, BV4) afin de respecter les contraintes de pureté tout en maximisant les performances.
* **Analyse comparative :** J'ai calculé les indicateurs de performance clés (productivité, ratio W/F, facteurs de dilution) pour chaque débit. Cette analyse a permis de conclure que le **débit optimal pour la séparation Sel/Glucose est de 2 BV/h**, offrant le meilleur compromis entre une productivité élevée et une faible consommation d'éluant.

***

## 🚀 Utilisation des Scripts

1.  **Simulation d'Élution Simple (`Simu_Chromato_Batch_EI.m`)**
    * Ouvrez le script dans MATLAB.
    * Dans la section "Chargement des paramètres", modifiez les valeurs des variables (ex: `NETcol_AF`, `Porosite_AF`, `Ks_A...`) pour ajuster le modèle.
    * Assurez-vous que le fichier de données expérimentales (ex: `Exp_NaCl_100gL_1BVh`) est correctement nommé et placé dans le bon répertoire.
    * Exécutez le script pour visualiser la superposition des profils simulés et expérimentaux.

2.  **Simulation du Procédé ISMB (`Simu_ISMB_2F_4Z_EI.m`)**
    * Ouvrez le script dans MATLAB.
    * Injectez les paramètres de rétention et de dispersion estimés à l'étape précédente pour les composés à séparer.
    * Ajustez les volumes `BV1_PT`, `BV2_PT`, `BV3_PT`, et `BV4_PT` de manière itérative.
    * Exécutez la simulation (qui peut prendre plusieurs cycles pour atteindre le régime permanent) pour visualiser les profils de concentration dans les 4 zones et les performances de séparation.
