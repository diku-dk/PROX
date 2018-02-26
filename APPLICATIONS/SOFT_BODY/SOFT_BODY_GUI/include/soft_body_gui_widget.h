#ifndef SOFT_BODY_GUI_WIDGET_H
#define SOFT_BODY_GUI_WIDGET_H

#include <QTimer>
#include <QGLWidget>
#include <QApplication>
#include <QMouseEvent>

namespace soft_body
{
  namespace gui
  {

    class Widget
    : public QGLWidget
    {
      Q_OBJECT

    private:

      //QOpenGLContext m__context;  // QWindow way

      QTimer        m_timer;

    public:

      Widget( QGLFormat const & format, QWidget *parent = 0);
      ~Widget();

      QSize minimumSizeHint() const;
      QSize sizeHint() const;

      public slots:

      void animate();

      //signals:

    protected:

      void initializeGL();
      void paintGL();
      void resizeGL(int width, int height);
      void mousePressEvent(QMouseEvent *event);
      void mouseMoveEvent(QMouseEvent *event);
      void mouseReleaseEvent(QMouseEvent *event);
      void keyPressEvent(QKeyEvent *e);
      
    };
    
  }// end namespace gui
  
}// end namespace soft_body

// SOFT_BODY_GUI_WIDGET
#endif
