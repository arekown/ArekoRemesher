# User interface texts. English is the source language in the code;
# German comes from this table via bpy.app.translations.

import bpy

from .hilfe import HELP

# Blender looks up labels and panel titles in the default context "*"
# and operator names in "Operator"; texts used in both are listed for both.
_CONTEXTS_FOR_BUTTONS = ("*", "Operator")

_GERMAN = {
    "Core ready": "Kern bereit",
    "Core not available": "Kern nicht verfügbar",
    "Create a quad mesh from the selected object": "Erzeugt ein Quad-Netz aus dem ausgewählten Objekt",
    "Select exactly one object": "Genau ein Objekt auswählen",
    "The selected object is not a mesh": "Das ausgewählte Objekt ist kein Netz",
    "{} Quads": "{} Vierecke",
    # insight
    "Insight": "Einblick",
    "Shows a stage as its own object.": "Zeigt eine Stufe als eigenes Objekt.",
    "Shows a stage as its own object": "Zeigt eine Zwischenstufe als eigenes Objekt",
    "Removes the stage object and shows the original again": "Entfernt das Stufenobjekt und zeigt das Original wieder",
    "Stage: {}": "Stufe: {}",
    "{} Points": "{} Punkte",
    "{} Parts": "{} Teile",
    "Part {}: {} Quads": "Teil {}: {} Vierecke",
    "Part {}: skipped": "Teil {}: übersprungen",
    "Part Color {}": "Teilfarbe {}",
    "Size Color {}": "Größenfarbe {}",
    "Grid Color {}": "Gitterfarbe {}",
    "No Field": "Ohne Feld",
    "Quad Faces": "Vierecke",
    "Triangle Faces": "Dreiecke",
    "Faces with 5+ Corners": "Fünf- und Mehrecke",
    # edge flow
    "Edge Flow": "Kantenfluss",
    "Guide Curves": "Führungskurven",
    "Kind for New Curves": "Art für neue Kurven",
    "Adds the selected curve objects as guide curves": "Fügt die ausgewählten Kurvenobjekte als Führungskurven hinzu",
    "Removes the chosen guide curve from the list": "Entfernt die gewählte Führungskurve aus der Liste",
    "Quads beside the curve follow its direction": "Vierecke neben der Kurve folgen ihrer Richtung",
    "An edge loop of the result follows the curve": "Eine Kantenschleife des Ergebnisses folgt der Kurve",
    "An edge loop follows the curve and stays sharp": "Eine Kantenschleife folgt der Kurve und bleibt scharf",
    "The curve is ignored": "Die Kurve wird übergangen",
    # symmetry
    "Symmetry": "Symmetrie",
    "Mirror on": "Spiegeln an",
    "Computed Half": "Gerechnete Hälfte",
    "Mirrored Half": "Gespiegelte Hälfte",
    # messages of the core status codes
    "Invalid input": "Ungültige Eingabe",
    "Out of memory": "Zu wenig Speicher",
    "Internal error": "Interner Fehler",
    "No mesh": "Kein Netz",
    "No result": "Kein Ergebnis",
    "Cancelled": "Abgebrochen",
    "Invalid settings": "Ungültige Einstellungen",
    "Stage not available yet": "Stufe noch nicht verfügbar",
    # create quads, run and messages
    "Done: {} quads": "Fertig: {} Vierecke",
    "{} Quads · {} Triangles": "{} Vierecke · {} Dreiecke",
    "Time: {} s": "Dauer: {} s",
    "Flat shading with Follow Normal Splits usually gives a poor result":
        "Flache Schattierung mit „Normalen-Knicken folgen“ ergibt meist ein schlechtes Ergebnis",
    "Flat shading and": "Flache Schattierung und",
    "Follow Normal Splits:": "Normalen-Knicken folgen:",
    "usually a poor result": "meist schlechtes Ergebnis",
    "Stops the running computation; nothing is created": "Bricht die laufende Rechnung ab; es entsteht nichts",
    # density
    "0 % = the same size everywhere": "0 % = überall gleich groß",
    "4 = finer · 0.25 = coarser": "4 = feiner · 0,25 = gröber",
    "Vertex paint with the colour of Brush Density: red gives smaller quads, cyan larger ones":
        "Vertex-Malen mit der Farbe der Pinsel-Dichte: Rot ergibt kleinere Vierecke, Cyan größere",
    "Density": "Dichte",
    # face groups, reset, shortcut
    "Gives every face set of the selected object a new material in a colour of its own":
        "Gibt jeder Flächengruppe des ausgewählten Objekts ein neues Material in eigener Farbe",
    "Face Group {}": "Flächengruppe {}",
    "No face groups found": "Keine Flächengruppen gefunden",
    "{} materials created": "{} Materialien angelegt",
    "Sets all settings back to their defaults": "Setzt alle Einstellungen auf ihre Vorgaben zurück",
    "Settings reset": "Einstellungen zurückgesetzt",
    "Shortcut: {}": "Tastenkürzel: {}",
}
# Panel titles live in their own context: Blender's own catalogue would win over
# ours for some of them.
PANEL_CONTEXT = "areko_panel"
_GERMAN_PANELS = {
    "Create Quads": "Quads erzeugen",
    "Density": "Dichte",
    "Fine Control": "Feinsteuerung",
    "Edge Flow": "Kantenfluss",
    "Guide Curves": "Führungskurven",
    "Borders": "Ränder",
    "Seams": "Nähte",
    "Symmetry": "Symmetrie",
    "Output": "Ausgabe",
    "Insight": "Einblick",
}
# Tooltips of the settings (property descriptions), registered in the default
# context and in the settings context.
_GERMAN_TOOLTIPS = {
    "About how many quads the result gets": "Ungefähr so viele Vierecke bekommt das Ergebnis",
    "Number of quads as a share of the input's faces": "Anzahl der Vierecke als Anteil an den Flächen der Eingabe",
    "Length of the quad edges, in object units": "Länge der Viereckskanten, in Objekteinheiten",
    "Off: the number of quads may change to suit the shape. On: it is kept":
        "Aus: Die Anzahl darf sich der Form anpassen. An: Sie wird eingehalten",
    "The first colour attribute sets the size: red gives smaller quads, cyan larger ones":
        "Das erste Farbattribut bestimmt die Größe: Rot ergibt kleinere Vierecke, Cyan größere",
    "The point attribute \"density\" sets the size: 4 gives smaller quads, 0.25 larger ones":
        "Das Punktattribut „density“ bestimmt die Größe: 4 ergibt kleinere Vierecke, 0,25 größere",
    "How much smaller the quads get where the surface bends; 0 = the same size everywhere":
        "Wie viel kleiner die Vierecke werden, wo sich die Fläche biegt; 0 = überall gleich groß",
    "Quads around a full bend; 0 = from Curvature Densifies":
        "Vierecke um eine volle Rundung; 0 = aus „Krümmung verdichtet“",
    "Ratio of the largest to the smallest quads; 0 = from Curvature Densifies, 1 = one size":
        "Verhältnis der größten zu den kleinsten Vierecken; 0 = aus „Krümmung verdichtet“, 1 = eine Größe",
    "Largest ratio of the long to the short side of a quad":
        "Größtes Verhältnis der langen zur kurzen Seite eines Vierecks",
    "Edge loops follow the borders between materials": "Kantenschleifen folgen den Grenzen zwischen Materialien",
    "Edge loops follow the borders between face sets": "Kantenschleifen folgen den Grenzen zwischen Flächengruppen",
    "Edge loops follow the borders between smoothing groups":
        "Kantenschleifen folgen den Grenzen zwischen Glättungsgruppen",
    "Edge loops follow the edges marked sharp": "Kantenschleifen folgen den als scharf markierten Kanten",
    "Edge loops follow the edges where the normals split":
        "Kantenschleifen folgen den Kanten, an denen die Normalen knicken",
    "Finds creases where the faces meet at a sharp angle": "Findet Knicke, wo Flächen in spitzem Winkel aufeinandertreffen",
    "Merges narrow strips between face groups, as CAD meshes have them":
        "Fasst schmale Streifen zwischen Flächengruppen zusammen, wie sie CAD-Netze haben",
    "Welds open seams whose border points lie exactly on top of each other":
        "Verschweißt offene Nähte, deren Randpunkte genau aufeinanderliegen",
    "Welded seams count as creases": "Verschweißte Nähte zählen als Knick",
    "Smooths out noise of the surface before the computation": "Glättet Rauschen der Fläche vor der Rechnung",
    "Faces meeting at a sharper angle than this form a crease":
        "Flächen, die spitzer als dieser Winkel aufeinandertreffen, bilden einen Knick",
    "Edge loops follow the open borders": "Kantenschleifen folgen den offenen Rändern",
    "The points of the open borders stay where they are": "Die Punkte der offenen Ränder bleiben, wo sie sind",
    "Scales the size of the quads along frozen borders; 1 = as elsewhere":
        "Skaliert die Größe der Vierecke an festgehaltenen Rändern; 1 = wie sonst",
    "The result is mirrored across the X plane": "Das Ergebnis wird an der X-Ebene gespiegelt",
    "The result is mirrored across the Y plane": "Das Ergebnis wird an der Y-Ebene gespiegelt",
    "The result is mirrored across the Z plane": "Das Ergebnis wird an der Z-Ebene gespiegelt",
    "The mirror plane runs through the middle of the object; off: through its origin":
        "Die Spiegelebene geht durch die Mitte des Objekts; aus: durch seinen Ursprung",
    "Finds the symmetry of every part and computes one half":
        "Findet die Symmetrie jedes Teils und rechnet eine Hälfte",
    "Holds the border at the mirror plane in place; only useful with Freeze Borders":
        "Hält den Rand an der Spiegelebene fest; nur mit „Ränder festhalten“ sinnvoll",
    "Splits the points of the result along its creases": "Trennt die Punkte des Ergebnisses entlang seiner Knicke",
    "The creases of the result get full edge crease": "Die Knicke des Ergebnisses bekommen vollen Kanten-Crease",
    "Dissolves leftover edges and separates points where the result is not manifold":
        "Löst übrige Kanten auf und trennt Punkte, wo das Ergebnis nicht mannigfaltig ist",
    "The colour attribute of the input is carried over to the result":
        "Das Farbattribut der Eingabe geht auf das Ergebnis über",
    "Cores that compute the separate pieces of the mesh side by side; 0 = all. The result is the same for every number":
        "Rechenkerne, die getrennte Teile des Netzes nebeneinander rechnen; 0 = alle. Das Ergebnis ist für jede "
        "Anzahl gleich",
    "Hides the input object once the result is there": "Blendet das Eingabeobjekt aus, sobald das Ergebnis da ist",
    "Keeps the prepared mesh for the next run: if only settings that act after the preparation change, it is "
    "skipped. The result is the same":
        "Behält das vorbereitete Netz für den nächsten Lauf: ändern sich nur Regler, die erst nach der Vorbereitung "
        "wirken, entfällt sie. Das Ergebnis ist gleich",
    "What sets the size of the quads": "Was die Größe der Vierecke bestimmt",
    "The result gets about this many quads": "Das Ergebnis bekommt ungefähr so viele Vierecke",
    "Colour of the density brush: 4 (red) gives smaller quads, 0.25 (cyan) larger ones":
        "Farbe des Dichte-Pinsels: 4 (rot) ergibt kleinere Vierecke, 0,25 (cyan) größere",
    # guide curves
    "The curve object that guides the edge flow": "Das Kurvenobjekt, das den Kantenfluss führt",
    "How the curve guides the result": "Wie die Kurve das Ergebnis führt",
    "On: the curve guides the result. Off: it is ignored":
        "An: Die Kurve führt das Ergebnis. Aus: Sie wird übergangen",
    "The width follows from the length of the curve": "Die Breite ergibt sich aus der Länge der Kurve",
    "How far beside the curve the quads follow its direction, in object units":
        "Wie weit neben der Kurve die Vierecke ihrer Richtung folgen, in Objekteinheiten",
    "On: the curves of the list guide the result. Off: they are ignored":
        "An: Die Kurven der Liste führen das Ergebnis. Aus: Sie werden übergangen",
    "Kind that added curves get": "Art, die hinzugefügte Kurven bekommen",
}
# Names of the settings (scene properties) live in their own context: Blender's
# own catalogue would win over ours for some of them (e.g. "Crease Angle").
SETTINGS_CONTEXT = "areko_settings"
_GERMAN_SETTINGS = {
    "Quads (Count)": "Vierecke (Anzahl)",
    "Quads (Share of Input, %)": "Vierecke (Anteil der Eingabe, %)",
    "Edge Length": "Kantenlänge",
    "Keep Count Strictly": "Anzahl streng einhalten",
    "Density from Color": "Dichte aus Farbe",
    "Density from Attribute": "Dichte aus Attribut",
    "Curvature Densifies": "Krümmung verdichtet",
    "Points per Circle": "Punkte je Rundung",
    "Density Range": "Dichte-Spanne",
    "Largest Side Ratio": "Größtes Seitenverhältnis",
    "Follow Material Borders": "Materialgrenzen folgen",
    "Follow Face Groups": "Flächengruppen folgen",
    "Follow Smoothing Groups": "Glättungsgruppen folgen",
    "Follow Sharp Edges": "Scharfen Kanten folgen",
    "Follow Normal Splits": "Normalen-Knicken folgen",
    "Detect Creases": "Knicke selbst erkennen",
    "Merge CAD Faces": "CAD-Flächen zusammenfassen",
    "Weld Open Seams": "Offene Nähte verschweißen",
    "Welded Seams as Crease": "Schweißnähte als Knick",
    "Denoise": "Entrauschen",
    "Crease Angle": "Knickwinkel",
    "Follow Borders": "Rändern folgen",
    "Freeze Borders": "Ränder festhalten",
    "Border Density Factor": "Randdichte-Faktor",
    "Mirror X": "Spiegeln X",
    "Mirror Y": "Spiegeln Y",
    "Mirror Z": "Spiegeln Z",
    "Mirror at Object Center": "Spiegeln an der Objektmitte",
    "Find Symmetry": "Symmetrie selbst finden",
    "Keep Half": "Halbseite festhalten",
    "Split at Creases": "An Knicken auftrennen",
    "Creases as Edge Crease": "Knicke als Kanten-Crease",
    "Clean Output": "Ausgabe säubern",
    "Keep Colors": "Farben übernehmen",
    "Cores": "Rechenkerne",
    "Hide Input": "Eingabe ausblenden",
    "Remember Preparation": "Vorbereitung merken",
    "Size by": "Größe nach",
    "Count": "Anzahl",
    "Share (%)": "Anteil (%)",
    "Brush Density": "Pinsel-Dichte",
    # guide curves
    "Curve": "Kurve",
    "Kind": "Art",
    "On": "An",
    "Automatic Width": "Breite automatisch",
    "Width": "Breite",
    "Use Guide Curves": "Führungskurven benutzen",
    "Kind for New Curves": "Art für neue Kurven",
    "Band": "Band",
    "Edge Loop": "Kantenschleife",
    "Crease": "Knick",
    "Off": "Aus",
    "Quads beside the curve follow its direction": "Vierecke neben der Kurve folgen ihrer Richtung",
    "An edge loop of the result follows the curve": "Eine Kantenschleife des Ergebnisses folgt der Kurve",
    "An edge loop follows the curve and stays sharp": "Eine Kantenschleife folgt der Kurve und bleibt scharf",
    "The curve is ignored": "Die Kurve wird übergangen",
    "One entry of the guide curve list": "Ein Eintrag der Führungskurven-Liste",
    "Guide curves of a scene": "Führungskurven einer Szene",
    "List of the guide curves": "Liste der Führungskurven",
    "One part of the shown stage.": "Ein Teil der gezeigten Stufe.",
}
_GERMAN_BUTTONS = {
    "Create Quads": "Quads erzeugen",
    "Show Stage": "Stufe zeigen",
    "Close Stage": "Stufe schließen",
    "Stage": "Stufe",
    "Add Curves": "Kurven hinzufügen",
    "Remove Curve": "Kurve entfernen",
    "Cancel (Esc)": "Abbrechen (Esc)",
    "Reset": "Zurücksetzen",
    "Face Groups to Materials": "Flächengruppen zu Materialien",
    "Paint Density": "Dichte malen",
    "Help": "Hilfe",
}
# Stage names live in their own context: several of them (e.g. "Input") already
# have a different meaning in Blender's own catalogue.
STAGE_CONTEXT = "areko_stage"
_GERMAN_STAGES = {
    "Stage": "Stufe",
    "Input": "Eingang",
    "Prepared": "Vorbereitet",
    "Cleaned": "Bereinigt",
    "Mirror Half": "Spiegelhälfte",
    "Feature Lines": "Merkmalslinien",
    "Parts": "Teile",
    "Curves": "Kurven",
    "Curvature": "Krümmung",
    "Cross Field": "Kreuzfeld",
    "Cuts": "Schnitte",
    "Size Field": "Größenfeld",
    "Parameterization": "Parametrisierung",
    "Raw Mesh": "Rohnetz",
    "Quads": "Vierecke",
    "Result": "Ergebnis",
}
# Help behind the "?" buttons, per topic of hilfe.HELP (the English text is the key there).
_GERMAN_HELP = {
    "size": "Legt fest, wie fein das Ergebnis wird. „Anzahl“ zielt auf ungefähr so viele Vierecke, „Anteil“ nimmt "
            "einen Prozentsatz der Flächen der Eingabe, „Kantenlänge“ gibt die Länge der Viereckskanten in "
            "Objekteinheiten vor. Nimm die Anzahl für ein festes Budget, etwa bei einem Spielobjekt, und die "
            "Kantenlänge, wenn mehrere Objekte gleich große Vierecke bekommen sollen.",
    "exact_count": "Aus: Die Zahl der Vierecke darf sich etwas der Form anpassen. An: Sie wird streng eingehalten. "
                   "Schalte es ein, wenn ein festes Budget wichtiger ist als der Kantenfluss, etwa für eine "
                   "vereinfachte Fassung mit fester Obergrenze.",
    "create_quads": "Berechnet aus dem ausgewählten Netzobjekt samt seinen Modifikatoren ein neues Vierecksnetz; "
                    "das ursprüngliche Netz bleibt unverändert. Blender bleibt währenddessen bedienbar, Esc bricht "
                    "ab. Das Ergebnis ist ein neues Objekt mit dem Namen der Eingabe. Tastenkürzel im "
                    "3D-Ansichtsfenster: Strg Umschalt Q.",
    "cancel_quads": "Bricht die laufende Rechnung sofort ab. Es entsteht nichts, und die Eingabe bleibt, wie sie "
                    "war. Nutze es, wenn ein Lauf zu lange dauert oder du vorher noch etwas einstellen willst.",
    "curvature": "Macht die Vierecke kleiner, wo sich die Fläche stark biegt, und größer auf flachen Bereichen. "
                 "Erhöhe den Wert bei Figuren und organischen Formen, damit Ohren, Finger oder Falten genug "
                 "Vierecke bekommen. 0 ergibt überall etwa gleich große Vierecke, etwa für eine Stoffsimulation.",
    "attribute_density": "Liest das Punktattribut „density“ (Kommazahl) der Eingabe: Werte bis 4 ergeben kleinere "
                         "Vierecke, Werte bis 0,25 größere, 1 ändert nichts. Nutze es, wenn du die Dichte selbst "
                         "berechnest, etwa mit Geometry Nodes. Von Hand geht es mit „Dichte malen“ leichter.",
    "color_density": "Liest das erste Farbattribut der Eingabe als Dichtekarte: Rot ergibt kleinere Vierecke, Cyan "
                     "größere, Weiß ändert nichts. Schalte es ein, um einzelnen Bereichen mehr oder weniger "
                     "Vierecke zu geben, etwa dem Gesicht einer Figur. „Dichte malen“ schaltet es selbst ein.",
    "density_paint": "Farbe des Pinsels für „Dichte malen“: 4 (Rot) ergibt die kleinsten Vierecke, 0,25 (Cyan) die "
                     "größten, 1 (Weiß) ändert nichts. Du kannst den Wert beim Malen ändern; der Pinsel übernimmt "
                     "die neue Farbe sofort.",
    "paint_density": "Schaltet das aktive Objekt in die Punktmalerei mit der Farbe der Pinsel-Dichte und schaltet "
                     "„Dichte aus Farbe“ ein. Ein Netz ohne Farben bekommt zuerst ein weißes Farbattribut; "
                     "sonst wird in sein erstes Farbattribut gemalt, und die Farben darin werden überschrieben. Male "
                     "Rot, wo du mehr Details willst, und Cyan, wo weniger Vierecke reichen.",
    "points_per_circle": "Wie viele Vierecke um eine volle Rundung laufen, etwa einmal um einen Zylinder. 0 "
                         "übernimmt den Wert aus „Krümmung verdichtet“. Erhöhe ihn, wenn runde Teile wie Rohre oder "
                         "Finger zu eckig wirken.",
    "density_ratio": "Begrenzt, wie viel größer die größten Vierecke als die kleinsten sein dürfen. 0 übernimmt den "
                     "Wert aus „Krümmung verdichtet“, 1 macht alle Vierecke gleich groß. Erhöhe ihn, wenn große "
                     "flache Bereiche grob bleiben und kleine Details trotzdem feine Vierecke bekommen sollen.",
    "max_quad_ratio": "Wie lang ein Viereck im Verhältnis zu seiner Breite werden darf. Kleinere Werte halten die "
                      "Vierecke nahezu quadratisch, größere erlauben längere Vierecke, wo die Form es verlangt. Die "
                      "Vorgabe passt für die meisten Modelle.",
    "detect_creases": "Findet harte Kanten am Winkel zwischen benachbarten Flächen („Knickwinkel“) und lässt "
                      "Kantenschleifen ihnen folgen. Lass es bei technischen Modellen und CAD-Teilen an. Schalte es "
                      "bei glatten organischen Formen und Scans aus, wo kleine Knicke keine Schleifen anziehen "
                      "sollen.",
    "crease_angle": "Mit „Knicke selbst erkennen“ gilt als Knick, wo Flächen spitzer als dieser Winkel "
                    "aufeinandertreffen. Senke ihn, damit auch weichere Kanten zählen, etwa an abgeschrägten Teilen; "
                    "erhöhe ihn, wenn Schleifen auf unwichtigen kleinen Knicken landen.",
    "follow_materials": "Kantenschleifen laufen an den Grenzen zwischen Materialien entlang, sodass jedes Material "
                        "saubere eigene Kanten bekommt. Schalte es ein, wenn die Materialien echte Teile markieren, "
                        "etwa die Platten eines technischen Modells. Es wirkt nur, wenn das Objekt mindestens zwei "
                        "verschiedene Materialien hat.",
    "follow_face_groups": "Kantenschleifen laufen an den Grenzen zwischen Flächengruppen (Face Sets) entlang, "
                          "die du im Skulpturmodus zeichnest. Schalte es ein, um Schleifen genau dorthin zu legen, wo du sie "
                          "haben willst, etwa um ein Auge oder entlang einer Fuge, ohne die Materialien anzufassen.",
    "follow_smoothing_groups": "Kantenschleifen laufen an den Grenzen zwischen Glättungsgruppen entlang. Blender "
                               "selbst hat keine; gelesen wird das Flächenattribut „smoothing_groups“, das die "
                               "Ergebnisse dieser Erweiterung tragen. Schalte es ein, wenn du ein solches Ergebnis "
                               "noch einmal umrechnest und seine Grenzen behalten willst.",
    "follow_sharp_edges": "Kantenschleifen laufen an den Kanten entlang, die im Bearbeitungsmodus als scharf "
                          "markiert sind. Das ist der direkteste Weg, eine Schleife auf eine harte Kante zu legen, "
                          "etwa bei technischen Modellen und CAD-Teilen. Lass es aus, wenn die Markierungen nur der "
                          "Schattierung dienen.",
    "follow_normal_splits": "Kantenschleifen laufen an den Kanten entlang, an denen die Normalen knicken, etwa bei "
                            "eigenen Normalen importierter Spiel- oder CAD-Modelle. Nutze es, wenn die harten "
                            "Kanten nur in den Normalen stecken. Bei flacher Schattierung knickt jede Kante, das "
                            "ergibt meist ein schlechtes Ergebnis.",
    "denoise": "Glättet kleines Rauschen der Fläche, bevor die Vierecke berechnet werden; scharfe Kanten bleiben. "
               "Lass es bei Scans und im Skulpturmodus geformten Modellen an, dort wird der Fluss ruhiger. Schalte "
               "es aus, wenn feine Details eines sauberen Modells den Fluss lenken sollen.",
    "face_groups_to_materials": "Gibt jeder Flächengruppe des ausgewählten Objekts ein neues Material in eigener "
                                "Farbe; jede Fläche nutzt danach das Material ihrer Gruppe. So siehst du die "
                                "Flächengruppen auch außerhalb des Skulpturmodus, und „Materialgrenzen folgen“ kann "
                                "sie nutzen. Es braucht mindestens zwei Flächengruppen.",
    "use_guide_curves": "Schaltet alle Führungskurven der Liste auf einmal an oder aus. Schalte es aus, um das "
                        "Ergebnis mit und ohne Kurven zu vergleichen, ohne die Liste zu verlieren.",
    "guide_curve_list": "Kurvenobjekte, die den Kantenfluss lenken. Wähle Kurvenobjekte aus und drücke +, um sie "
                        "hinzuzufügen; - entfernt den gewählten Eintrag, und das Häkchen eines Eintrags schaltet "
                        "diese Kurve an oder aus. Zeichne die Kurven auf oder dicht über der Fläche, etwa um einen "
                        "Hals oder entlang einer Fuge.",
    "curve": "Das Kurvenobjekt des gewählten Eintrags. Wähle hier ein anderes, um die Kurve zu tauschen und Art und "
             "Breite zu behalten.",
    "kind": "Wie die Kurve das Ergebnis lenkt: Bei „Band“ folgen die Vierecke neben ihr ihrer Richtung, "
            "„Kantenschleife“ legt eine Kantenschleife auf sie, „Knick“ ebenso und hält diese Schleife scharf, „Aus“ "
            "übergeht sie. Nimm „Kantenschleife“ für genaue Schleifen, etwa um einen Mund, und „Band“ für eine grobe "
            "Richtung, etwa entlang eines Arms.",
    "automatic_width": "An: Die Breite, in der die Kurve wirkt, ergibt sich aus ihrer Länge. Schalte es aus, um die "
                       "Breite selbst einzustellen, etwa wenn ein Band einen breiteren oder schmaleren Streifen "
                       "lenken soll.",
    "width": "Wie weit neben der Kurve die Vierecke ihrer Richtung folgen, in Objekteinheiten. Es zählt nur, wenn "
             "„Breite automatisch“ aus ist. Größer gibt einen breiten, ruhigen Fluss; kleiner, wenn die Kurve nur "
             "dicht an sich selbst wirken soll.",
    "new_kind": "Die Art, die Kurven beim Hinzufügen mit + bekommen. Stelle sie ein, bevor du mehrere Kurven "
                "derselben Art hinzufügst; jeder Eintrag lässt sich danach noch ändern.",
    "follow_borders": "Kantenschleifen laufen an den offenen Rändern des Netzes entlang, etwa an Ärmelöffnungen, "
                      "Säumen oder dem Rand eines Lochs. Lass es bei Kleidung und anderen offenen Formen an; schalte "
                      "es nur aus, wenn die Form des Randes den Fluss nicht lenken soll.",
    "freeze_borders": "Die Punkte der offenen Ränder bleiben genau, wo sie sind, und die Größe der Vierecke richtet "
                      "sich dann nach dem Abstand dieser Punkte, nicht nach der Anzahl. Nutze es, wenn das Ergebnis "
                      "am Rand genau an ein anderes Netz passen muss, etwa ein Kleidungsteil, das an einen Körper "
                      "anschließt.",
    "border_density_factor": "Nur mit „Ränder festhalten“: skaliert die Größe der Vierecke, die sich aus dem "
                             "Abstand der festgehaltenen Randpunkte ergibt. 1 behält diesen Abstand, Werte über 1 "
                             "ergeben größere Vierecke, Werte unter 1 kleinere.",
    "weld_seams": "Verbindet offene Nähte, deren Randpunkte genau aufeinanderliegen, sodass das Netz dort "
                  "geschlossen ist. Schalte es ein bei Modellen, die geschlossen aussehen, aber an Nähten getrennt "
                  "sind, etwa bei manchen Importen; sonst zählt jede Naht als offener Rand.",
    "weld_as_crease": "Mit „Offene Nähte verschweißen“ verbundene Nähte zählen als Knick: Kantenschleifen folgen "
                      "ihnen, und sie bleiben scharf. Nutze es, wenn die Nähte echte Kanten sind, etwa wo die "
                      "Platten eines CAD-Modells zusammentreffen; lass es aus, wenn sie nur eine technische Trennung "
                      "waren.",
    "merge_cad_groups": "Schlägt schmale Streifen zwischen Flächengruppen ihren Nachbarn zu, etwa die kleinen "
                        "Rundungen aus CAD-Exporten, damit sie keine eigenen Kantenschleifen anziehen. Es wirkt auf "
                        "die Grenzen, denen der Kantenfluss folgt (Materialien, Flächengruppen, Glättungsgruppen).",
    "mirror_axes": "Berechnet eine Hälfte und spiegelt sie an den gewählten Ebenen, sodass das Ergebnis genau "
                   "symmetrisch ist. X passt für eine Figur, die nach vorn schaut; mehrere Achsen lassen sich "
                   "kombinieren. Wo die Ebenen liegen, bestimmt „Spiegeln an der Objektmitte“.",
    "mirror_at_center": "An: Die Spiegelebene geht durch die Mitte des Objekts, gemessen an seinen äußeren Maßen; "
                        "aus: durch seinen Ursprung. Schalte es aus, wenn das Objekt um seinen Ursprung symmetrisch "
                        "ist, sein Umriss aber nicht, etwa eine Figur, die etwas in einer Hand hält.",
    "detect_symmetry": "Findet die Spiegelebene jedes Teils selbst und berechnet von jedem eine Hälfte, ohne "
                       "Achsenangabe. Nutze es für Teile, die symmetrisch, aber nicht an den Achsen ausgerichtet "
                       "sind, oder für mehrere Teile mit eigenen Ebenen. Es vergleicht die Flächen beider Seiten und "
                       "klappt deshalb am besten bei modellierten Netzen.",
    "freeze_half": "Für ein Netz, das nur eine Hälfte ist, etwa für einen Modifikator „Spiegeln“: Der offene Rand an "
                   "der Spiegelebene bleibt, wo er ist, damit beide Hälften weiter zusammenpassen, und das Ergebnis "
                   "bleibt eine Hälfte. Wähle die Achse unter „Spiegeln an“; sinnvoll nur zusammen mit „Ränder "
                   "festhalten“.",
    "hide_input": "Blendet das Eingabeobjekt aus, sobald das Ergebnis da ist, damit die beiden nie übereinanderliegen. "
                  "Schalte es aus, um beide sichtbar zu lassen; schieb dann eines zur Seite, um sie zu vergleichen.",
    "split_at_creases": "Trennt die Punkte des Ergebnisses entlang seiner Knicke, sodass die Flächen auf beiden "
                        "Seiten sie nicht mehr teilen. Nutze es für Exporte, bei denen harte Kanten echte Trennungen "
                        "sein müssen, etwa zum Backen von Texturen oder für eine Spiel-Engine; lass es aus, wenn du "
                        "am Ergebnis weitermodellierst.",
    "creases_to_edge_crease": "Gibt den Knicken des Ergebnisses vollen Kanten-Crease, sodass ein Modifikator "
                              "„Oberflächenunterteilung“ sie scharf hält. Schalte es ein bei technischen Modellen, "
                              "die du danach unterteilst.",
    "clean_output": "Räumt das Ergebnis auf: Es löst übrige Kanten auf und trennt Punkte, wo das Netz nicht sauber "
                    "zusammenhängt (nicht mannigfaltig). Schalte es ein, wenn das Ergebnis in Programme geht, die "
                    "ein sauberes Netz brauchen, etwa eine Simulation oder den 3D-Druck.",
    "keep_colors": "Überträgt das Farbattribut der Eingabe auf das Ergebnis. Schalte es ein, wenn die Farben später "
                   "gebraucht werden, etwa gemalte Punktfarben; lass es aus, wenn sie nur als Dichtekarte dienten.",
    "thread_count": "Wie viele Rechenkerne die getrennten Teile des Netzes nebeneinander berechnen; 0 nutzt alle. "
                    "Es beschleunigt nur Netze aus mehreren Teilen, und das Ergebnis ist bei jeder Anzahl gleich. "
                    "Senke den Wert, wenn andere Programme flüssig bleiben sollen.",
    "remember_preparation": "Behält das vorbereitete Netz nach einem Lauf. Änderst du danach nur Einstellungen, die "
                            "später wirken, etwa „Krümmung verdichtet“ oder die Ausgabe, überspringt der nächste "
                            "Lauf die Vorbereitung und ist schneller; das Ergebnis ist gleich. Schalte es aus, um "
                            "bei sehr großen Netzen den Speicher freizugeben.",
    "reset": "Setzt alle Einstellungen des Panels auf ihre Vorgaben zurück, auch „Führungskurven benutzen“ und „Art "
             "für neue Kurven“; die Liste der Führungskurven bleibt. Nutze es für einen sauberen Neuanfang, wenn "
             "eine Mischung von Einstellungen seltsame Ergebnisse gibt.",
    "show_stage": "Rechnet bis zur gewählten Stufe und zeigt sie als eigenes Objekt an der Stelle der Eingabe, die "
                  "so lange ausgeblendet ist. So siehst du, warum ein Ergebnis so aussieht: „Merkmalslinien“ zeigt "
                  "etwa die gefundenen Knicke und „Größenfeld“, wo die Vierecke klein werden. Geänderte "
                  "Einstellungen frischen die gezeigte Stufe nach einem Moment auf.",
    "close_stage": "Entfernt das Stufenobjekt und zeigt die Eingabe wieder. „Quads erzeugen“ macht das vor dem "
                   "Start von selbst.",
}


class Translations:
    """Registers the translation table under the add-on's name."""

    LOCALES = ("de_DE",)

    def __init__(self, owner: str):
        self._owner = owner

    @staticmethod
    def _table() -> dict:
        german = {("*", en): de for en, de in _GERMAN.items()}
        for en, de in _GERMAN_BUTTONS.items():
            for context in _CONTEXTS_FOR_BUTTONS:
                german[(context, en)] = de
        for en, de in _GERMAN_SETTINGS.items():
            german[(SETTINGS_CONTEXT, en)] = de
        for en, de in _GERMAN_PANELS.items():
            german[(PANEL_CONTEXT, en)] = de
        for en, de in _GERMAN_TOOLTIPS.items():
            german[("*", en)] = de
            german[(SETTINGS_CONTEXT, en)] = de
        for en, de in _GERMAN_STAGES.items():
            german[(STAGE_CONTEXT, en)] = de
        for topic, de in _GERMAN_HELP.items():
            german[("*", HELP[topic])] = de
        return {locale: german for locale in Translations.LOCALES}

    def register(self) -> None:
        bpy.app.translations.register(self._owner, self._table())

    def unregister(self) -> None:
        bpy.app.translations.unregister(self._owner)


def iface(text: str, context: str | None = None) -> str:
    """Translate a text built at runtime (labels with values)."""
    return bpy.app.translations.pgettext_iface(text, context)


def _german() -> bool:
    view = bpy.context.preferences.view
    return view.use_translate_interface and bpy.app.translations.locale.startswith("de")


def grouped(count: int) -> str:
    """Whole number with grouped digits: "75 000" in German, "75,000" in English."""
    return f"{count:,}".replace(",", " " if _german() else ",")


def decimal(value: float, digits: int) -> str:
    """Number with `digits` decimals: "3,4" in German, "3.4" in English."""
    text = f"{value:.{digits}f}"
    return text.replace(".", ",") if _german() else text


def report_text(text: str) -> str:
    """Translate an operator report."""
    return bpy.app.translations.pgettext_rpt(text)
