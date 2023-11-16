# Lines

Дисциплина: **Разработка САПР**

Тема: **Геометрический решатель**

```
pip install -r requirements.txt
```

Запуск
```
python main.py
```

Описание
```
class Constraint:
    name = 'point_belongs_line_constraint'
    objects: [{'type': 'line', 'obj': 0},{'type': 'point', 'obj': 1}]
    value = None

class Constraint:
    name = 'points_dist_constraint'
    objects: [{'type': 'point', 'obj': 0},{'type': 'point', 'obj': 1}]
    value = 50.
```

- словарик points у storage содержит точки в виде:
QPoint - https://doc.qt.io/qtforpython/PySide2/QtCore/QPoint.html
- есть два метода: x() и y()
- словарик line: ключ - id отрезок {'p1_id': first_id, 'p2_id': second_id}