class MonObjet:
    def __init__(self, value):
        self.ma_valeur = value


def fct(obj):
    return obj.ma_valeur == 2


ma_liste_objets = []
for i in range(10):
    ma_liste_objets.append(MonObjet(i))

print(ma_liste_objets)

mon_objet = list(filter(lambda obj: obj.ma_valeur == 6, ma_liste_objets))[0]
print(mon_objet.ma_valeur)

v = []
for p in ma_liste_objets:
    v.append(p.ma_valeur)
print(v)