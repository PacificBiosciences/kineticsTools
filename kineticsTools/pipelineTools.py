def consumer(func):
    def start(*args, **kwargs):
        c = func(*args, **kwargs)
        next(c)
        return c
    return start


def broadcast(source, consumers):
    for item in source:
        for c in consumers:
            c.send(item)


def concat(sources):
    for source in sources:
        for item in source:
            yield item
