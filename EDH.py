class Analitical:
    """
    Classe contendo as soluções analíticas.
    Suporta somente problemas 1D

    dimension :

    CI :

    CC : [[tipo da CC, tipo da CC],[lado esquerdo, lado direto]] -
        São as condições de contorno do problema. Podem ser do tipo Diritchlet (pressão) ou Newmman (vazão).

    grid

    """

    def __init__(self, dimension : int, CI:[], CC : [list, list], grid : tuple | list, ):
        self.dimension = dimension

class Numerical:
    def __init__(self):
        pass